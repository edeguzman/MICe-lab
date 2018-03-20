#ifdef _CH_
#pragma package <opencv>
#endif

#define CV_NO_BACKWARD_COMPATIBILITY

#ifndef _EiC
#include "cv.h"
#endif

#include "stdio.h"
#include "unistd.h"

//compile with: gcc -o CCimages -L/usr/local/lib -I/usr/local/include/opencv CCimages.c -lopencv_core -lopencv_highgui
//              gcc -o CCimages CCimages.c -L/usr/local/lib -I/usr/local/include/opencv -lopencv_core -lopencv_highgui `pkg-config --cflags --libs opencv`

int search_width=-1,search_height=-1;
int refx0,refy0,refx1,refy1;
int movx0,movy0,movx1,movy1;
int roix0,roiy0,roix1,roiy1;
int templ_width=-1,templ_height=-1;
int nCC_size=9;
float offset_guessx=0.0,offset_guessy=0.0;
double minVal,maxVal,meanVal;

char* refimgfile;
char* movimgfile;

int verbose=0,showimages=0,normCCresult=0;

IplImage *refimage = 0, *movimage=0, *subrefimage=0, *tempimage=0, *matchresultimg=0, *filtmatchresultimg=0;
IplImage *subrefimage32 = 0, *tempimage32 = 0;
CvPoint minLoc, maxLoc;
CvScalar max_local_avg;
int max_rect_size=5;

int Iclamp(int x, int a, int b)
{
    return x < a ? a : (x > b ? b : x);
}

IplImage* Sub_Image(IplImage *image, CvRect roi)
{
    IplImage *result;
    cvSetImageROI(image,roi);
    result = cvCreateImage( cvSize(roi.width, roi.height), image->depth, image->nChannels );
    cvCopy(image,result,NULL);
    cvResetImageROI(image); // release image ROI
    return result;
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
             search_width = atoi(optarg);
             break;
           case 'l':
             search_height = atoi(optarg);
             break;
           case 'x':
             offset_guessx = atof(optarg);
             break;
           case 'y':
             offset_guessy = atof(optarg);
             break;
           case 't':
             templ_width = atoi(optarg);
             break;
           case 'u':
             templ_height = atoi(optarg);
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

    if( (refimage = cvLoadImage( refimgfile, 1)) == 0 )
    {
        printf("Failed to load %s\n",refimgfile);
        return -1;
    }

    if( (movimage = cvLoadImage( movimgfile, 1)) == 0 )
    {
        printf("Failed to load %s\n",movimgfile);
        return -1;
    }

    //set search ROI
    if (search_width<0) search_width = (int) (0.8*((float) refimage->width));
    if (search_height<0) search_height = (int) (0.8*((float) refimage->height));
    if (templ_width<0) templ_width = 2+search_width/5;
    if (templ_height<0) templ_height = 2+search_height/5;

    if( (fabs(offset_guessx)>(float) refimage->width-templ_width)||
        (fabs(offset_guessy)>(float) refimage->height-templ_height) )
    {
        printf("Offset and/or template size incompatible with image dimensions.\n");
        return -1;
    }

    
    //set search area on reference image
    refx0 = refimage->width/2 + (int) offset_guessx - search_width/2 - templ_width/2;
    refx0 = Iclamp(refx0, 0, refimage->width-search_width-templ_width);
    refy0 = refimage->height/2 + (int) offset_guessy - search_height/2 - templ_height/2;
    refy0 = Iclamp(refy0, 0, refimage->height-search_height-templ_height);
    refx1 = search_width+templ_width;
    refy1 = search_height+templ_height;

    //set template area on image to place
    movx0 = movimage->width/2 - (int) offset_guessx - templ_width/2;
    movx0 = Iclamp(movx0, 0, movimage->width-templ_width);
    movy0 = movimage->height/2 - (int) offset_guessy - templ_height/2;
    movy0 = Iclamp(movy0, 0, movimage->height-templ_height);
    movx1 = templ_width;
    movy1 = templ_height;

    if (verbose){
        printf("templ_width=%d, templ_height=%d, search_width=%d, search_height=%d\n",templ_width,templ_height,search_width,search_height);
        printf("subrefimg: refx0=%d, refy0=%d, refx1=%d, refy1=%d\n",refx0,refy0,refx1,refy1);
        printf("tempimg  : movx0=%d, movy0=%d, movx1=%d, movy1=%d\n",movx0,movy0,movx1,movy1);
    }

    //extract a subimage from reference
    subrefimage = Sub_Image(refimage, cvRect( refx0, refy0, refx1, refy1 ) );
    cvSetImageCOI(subrefimage, 1); 
    cvMinMaxLoc(subrefimage,&minVal,&maxVal,&minLoc,&maxLoc,NULL);
    cvSetImageCOI(subrefimage, 0);
    subrefimage32 = cvCreateImage(cvGetSize(subrefimage), IPL_DEPTH_32F, subrefimage->nChannels);
    cvConvertScale(subrefimage, subrefimage32, 1.0/maxVal, 0.0);

    //extract a template for matching
    tempimage = Sub_Image(movimage, cvRect( movx0, movy0, movx1, movy1 ) );
    cvSetImageCOI(tempimage, 1); 
    cvMinMaxLoc(tempimage,&minVal,&maxVal,&minLoc,&maxLoc,NULL);
    cvSetImageCOI(tempimage, 0);
    tempimage32 = cvCreateImage(cvGetSize(tempimage), IPL_DEPTH_32F, tempimage->nChannels);
    cvConvertScale(tempimage, tempimage32, 1.0/maxVal, 0.0);

    //use cvMatchTemplate to find match
    //options CV_TM_SQDIFF,CV_TM_SQDIFF_NORMED,CV_TM_CCORR,CV_TM_CCORR_NORMED,CV_TM_CCOEFF,CV_TM_CCOEFF_NORMED
    matchresultimg = cvCreateImage(cvSize(subrefimage32->width - tempimage32->width + 1,\
                                          subrefimage32->height - tempimage32->height + 1),IPL_DEPTH_32F,1);
    cvMatchTemplate(subrefimage,tempimage,matchresultimg,CV_TM_CCOEFF_NORMED);
    //compare max to local average if -n
    if (normCCresult){
        filtmatchresultimg = cvCreateImage( cvSize(matchresultimg->width, matchresultimg->height), \
                                            matchresultimg->depth, matchresultimg->nChannels );
        cvSmooth(matchresultimg,filtmatchresultimg,CV_BLUR,nCC_size,nCC_size,0,0);
        cvAddWeighted(matchresultimg,(nCC_size*nCC_size+1.0)/(nCC_size*nCC_size),\
                      filtmatchresultimg,-1.0,0.0,matchresultimg);

        //cvAddWeighted(matchresultimg,(nCC_size*nCC_size-2.0)/(nCC_size*nCC_size-1.0),\
        //              filtmatchresultimg,-1.0*(nCC_size*nCC_size)/(nCC_size*nCC_size-1.0),0.0,filtmatchresultimg);
        //cvMul(matchresultimg,filtmatchresultimg,matchresultimg,1.0);
    }
    cvMinMaxLoc(matchresultimg,&minVal,&maxVal,&minLoc,&maxLoc,NULL);

    //mean intensity in images
    //meanVal = 0.5*( cvAvg(refimage,NULL) ).val[0] + 0.5*( cvAvg(movimage,NULL) ).val[0];
    meanVal = 0.5*( cvAvg(subrefimage,NULL) ).val[0] + 0.5*( cvAvg(tempimage,NULL) ).val[0];


    printf("Imean=%f, offset_x=%d, offset_y=%d, CC=%f\n",(float) meanVal,\
                                                         refx0-movx0+maxLoc.x+refimage->width/2-movimage->width/2,\
                                                         refy0-movy0+maxLoc.y+refimage->height/2-movimage->height/2,\
                                                         (float) maxVal );

    if (showimages)
    {
        char wndname[] = "Image";
        cvNamedWindow(wndname,0); 
        cvRectangle(refimage,cvPoint(refx0,refy0),\
                    cvPoint(refx0+search_width+templ_width,refy0+search_height+templ_height),\
                    cvScalar(255,255,255,0),1,8,0);
        cvShowImage(wndname,refimage);
        char wndname2[] = "Template";
        cvNamedWindow(wndname2,0);
        cvRectangle(movimage,cvPoint(movx0,movy0),\
                    cvPoint(movx0+templ_width,movy0+templ_height),\
                    cvScalar(255,255,255,0),1,8,0);
        cvShowImage(wndname2,movimage);
        char wndname3[] = "matchresult";
        cvNamedWindow(wndname3,0);
        cvNormalize( matchresultimg, matchresultimg, 1.0, 0.0, CV_MINMAX, NULL );
        cvShowImage(wndname3,matchresultimg);
        char wndname4[] = "SubImage";
        cvNamedWindow(wndname4,0); 
        cvShowImage(wndname4,subrefimage);
        char wndname5[] = "SubTemplate";
        cvNamedWindow(wndname5,0);
        cvShowImage(wndname5,tempimage);
        cvWaitKey(0);
    }
    return 0;
}
