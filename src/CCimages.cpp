#include <opencv2/opencv.hpp>
#include <getopt.h>

using namespace cv;
using namespace std;

int search_cols=-1,search_rows=-1;
int refx0,refy0,refx1,refy1;
int movx0,movy0,movx1,movy1;
int roix0,roiy0,roix1,roiy1;
int templ_cols=-1,templ_rows=-1;
int nCC_size=9;
float offset_guessx=0.0,offset_guessy=0.0;
double minVal,maxVal,meanVal;

char* refimgfile;
char* movimgfile;

int verbose=0,showimages=0,normCCresult=0;
int j;

Mat refimage, movimage, subrefimage, tempimage, matchresultimg, filtmatchresultimg;
vector<Mat> splitimg;

Point minLoc, maxLoc;
Scalar max_local_avg;
int max_rect_size=5;

int Iclamp(int x, int a, int b)
{
    return x < a ? a : (x > b ? b : x);
}

int main( int argc, char** argv )
{
    int c;
    char *filename;
    opterr = 0;
     
    while ((c = getopt (argc, argv, "x:y:w:l:t:u:vsn")) != -1)
         switch (c)
           {
           case 'w':
             search_cols = atoi(optarg);
             break;
           case 'l':
             search_rows = atoi(optarg);
             break;
           case 'x':
             offset_guessx = atof(optarg);
             break;
           case 'y':
             offset_guessy = atof(optarg);
             break;
           case 't':
             templ_cols = atoi(optarg);
             break;
           case 'u':
             templ_rows = atoi(optarg);
             break;
           case 'v':
             verbose = 1;
             break;
           case 's':
             showimages = 1;
             break;
           case 'n':
             normCCresult = 1;
             break;
           default:
             abort ();
           }

    if (argc-optind<2)
    {
        fprintf(stderr, "Two input images are required.\n");
        return 1;
    }
    refimgfile = argv[optind];
    movimgfile = argv[optind+1];

    refimage = imread( refimgfile, 1 );
    if( !refimage.data )
    {
        printf("Failed to load %s\n",refimgfile);
        return -1;
    }

    movimage = imread( movimgfile, 1 );
    if ( !movimage.data )
    {
        printf("Failed to load %s\n",movimgfile);
        return -1;
    }

    //set search ROI
    if (search_cols<0) search_cols = (int) (0.8*((float) refimage.cols));
    if (search_rows<0) search_rows = (int) (0.8*((float) refimage.rows));
    if (templ_cols<0) templ_cols = 2+search_cols/5;
    if (templ_rows<0) templ_rows = 2+search_rows/5;

    if( (fabs(offset_guessx)>(float) refimage.cols-templ_cols)||
        (fabs(offset_guessy)>(float) refimage.rows-templ_rows) )
    {
        printf("Offset and/or template size incompatible with image dimensions.\n");
        return -1;
    }

    
    //set search area on reference image
    refx0 = refimage.cols/2 + (int) offset_guessx - search_cols/2 - templ_cols/2;
    refx0 = Iclamp(refx0, 0, refimage.cols-search_cols-templ_cols);
    refy0 = refimage.rows/2 + (int) offset_guessy - search_rows/2 - templ_rows/2;
    refy0 = Iclamp(refy0, 0, refimage.rows-search_rows-templ_rows);
    refx1 = search_cols+templ_cols;
    refy1 = search_rows+templ_rows;

    //set template area on image to place
    movx0 = movimage.cols/2 - (int) offset_guessx - templ_cols/2;
    movx0 = Iclamp(movx0, 0, movimage.cols-templ_cols);
    movy0 = movimage.rows/2 - (int) offset_guessy - templ_rows/2;
    movy0 = Iclamp(movy0, 0, movimage.rows-templ_rows);
    movx1 = templ_cols;
    movy1 = templ_rows;

    if (verbose){
        printf("templ_cols=%d, templ_rows=%d, search_cols=%d, search_rows=%d\n",templ_cols,templ_rows,search_cols,search_rows);
        printf("subrefimg: refx0=%d, refy0=%d, refx1=%d, refy1=%d\n",refx0,refy0,refx1,refy1);
        printf("tempimg  : movx0=%d, movy0=%d, movx1=%d, movy1=%d\n",movx0,movy0,movx1,movy1);
    }

    //extract a subimage from reference
    subrefimage = refimage(Rect( refx0, refy0, refx1, refy1 ));

    //extract a template for matching
    tempimage = movimage(Rect( movx0, movy0, movx1, movy1 ));

    //use cvMatchTemplate to find match
    //options CV_TM_SQDIFF,CV_TM_SQDIFF_NORMED,CV_TM_CCORR,CV_TM_CCORR_NORMED,CV_TM_CCOEFF,CV_TM_CCOEFF_NORMED
    matchresultimg = Mat(subrefimage.rows - tempimage.rows+1, subrefimage.cols-tempimage.cols+1, CV_32FC1);
    matchTemplate(subrefimage,tempimage,matchresultimg,CV_TM_CCOEFF_NORMED);
    //compare max to local average if -n
    if (normCCresult){
        filtmatchresultimg = Mat( matchresultimg.size(), CV_32FC1 );
        GaussianBlur(matchresultimg,filtmatchresultimg,Size(nCC_size,nCC_size),0,0,BORDER_CONSTANT);
        addWeighted(matchresultimg,(nCC_size*nCC_size+1.0)/(nCC_size*nCC_size),filtmatchresultimg,-1.0,0.0,matchresultimg);
    }
    minMaxLoc(matchresultimg,&minVal,&maxVal,&minLoc,&maxLoc,Mat());

    //mean intensity in images
    meanVal = 0.5*( mean(subrefimage,Mat()) ).val[2] + 0.5*( mean(tempimage,Mat()) ).val[2];

    //output result to terminal
    printf("Imean=%f, offset_x=%d, offset_y=%d, CC=%f\n",(float) meanVal,\
                                                         refx0-movx0+maxLoc.x+refimage.cols/2-movimage.cols/2,\
                                                         refy0-movy0+maxLoc.y+refimage.rows/2-movimage.rows/2,\
                                                         (float) maxVal );

    if (showimages)
    {
        char wndname[] = "Image";
        namedWindow(wndname,0); 
        rectangle(refimage,Point(refx0,refy0),Point(refx0+search_cols+templ_cols,refy0+search_rows+templ_rows), \
                  Scalar(255,255,255,0),1,8,0);
        imshow(wndname,refimage);
        char wndname2[] = "Template";
        namedWindow(wndname2,0);
        rectangle(movimage,Point(movx0,movy0),Point(movx0+templ_cols,movy0+templ_rows), \
                  Scalar(255,255,255,0),1,8,0);
        imshow(wndname2,movimage);
        char wndname3[] = "matchresult";
        namedWindow(wndname3,0);
        normalize( matchresultimg, matchresultimg, 1.0, 0.0, NORM_MINMAX, -1, Mat() );
        imshow(wndname3,matchresultimg);
        char wndname4[] = "SubImage";
        namedWindow(wndname4,0); 
        imshow(wndname4,subrefimage);
        char wndname5[] = "SubTemplate";
        namedWindow(wndname5,0);
        imshow(wndname5,tempimage);
        cvWaitKey(0);
    }
    return 0;
}
