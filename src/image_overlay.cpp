#include <opencv2/opencv.hpp>
using namespace cv;

//compile with: 
// g++ image_overlay.cpp -o image_overlay -L/usr/local/lib -I/usr/local/include/opencv -L${exec_prefix}/lib -lopencv_shape -lopencv_stitching -lopencv_objdetect -lopencv_superres -lopencv_videostab -lopencv_calib3d -lopencv_features2d -lopencv_highgui -lopencv_videoio -lopencv_imgcodecs -lopencv_video -lopencv_photo -lopencv_ml -lopencv_imgproc -lopencv_flann -lopencv_viz -lopencv_core


void AddROI(Mat src, Mat overlay, Point location)
{
    int x,y,i;
    int ROIwidth,ROIheight;
    Mat srcroi,overlayroi;
    
    ROIwidth = (location.x+overlay.cols > src.cols) ? src.cols-location.x : overlay.cols;
    ROIheight = (location.y+overlay.rows > src.rows) ? src.rows-location.y : overlay.rows;
    if ((ROIwidth>0)&&(ROIheight>0))
    {
        srcroi = src(Rect( location.x, location.y, ROIwidth, ROIheight ));
        overlayroi = overlay(Rect( 0,0,ROIwidth,ROIheight));
        add(srcroi,overlayroi,srcroi);
    }
}

Mat DistMap(int width, int height, int nChannels)
{
    int x,y,midy,midx,norm;
    int ROIwidth,ROIheight;
    float cdist;
    Mat result;
    result = Mat( Size(width, height), ((nChannels==1) ? CV_32FC1 : CV_32FC3) );
    midx = width/2;
    midy = height/2;
    norm = (width>=height) ? width/2.0 : height/2.0;
    for (x=0; x<width; x++)
    {
        for (y=0; y<height; y++)
        {
            if (y<midy)
                cdist=(float) y; 
            else
                cdist=(float) (height-y-1); 
            if (x<midx)
                cdist=((float) x < cdist) ? (float) x : cdist; 
            else
                cdist=((float) (width-x-1) < cdist) ? (float) (width-x-1) : cdist; 
            if (nChannels==1)
                result.at<float>(Point(x,y)) = cdist/norm;
            else
                result.at<Vec3f>(Point(x,y)) = Vec<float,3>(cdist/norm,cdist/norm,cdist/norm);
        }
    }
    return result;
}

double ReturnTypeMax(int cvdepthtype)
{
    double datamax;
    switch (cvdepthtype)
    {
        case CV_32F:
        case CV_64F:
            datamax=1.0;     break;
        case CV_8U:
            datamax=255.0;   break;
        case CV_8S:
            datamax=128.0;   break;
        case CV_16U:
            datamax=65535.0; break;
        case CV_16S:
            datamax=32768.0;      break;
        case CV_32S:
            datamax=2147483648.0; break;
        default:
            datamax=1.0;
    }
    return datamax;
}

int output_width=0,output_height=0,j;
int outputtype=-1,outputdepth=-1;
int oddsizes=0;
int colformat=CV_LOAD_IMAGE_UNCHANGED;
int nChannels=1;
long int elementsize=4L;
long int imagebytesize;

char* nextimgfile;
char* overlayimgfile;

double minVal,maxVal;
double offset_x,offset_y;
double indatamax,outdatamax;
double outputscale=1.0;

Mat nextimage, overlayimage;
Mat next32, outimage;
Mat cdist, tdist;
Point offset;
Point minLoc, maxLoc;


int main( int argc, char** argv )
{
    int c;
    char *filename;
    opterr = 0;
     
    while ((c = getopt (argc, argv, "x:y:I:sigo")) != -1)
         switch (c)
           {
           case 'x':
             output_width = atoi(optarg);
             break;
           case 'y':
             output_height = atoi(optarg);
             break;
           case 'I':
             outputscale = atof(optarg);
             break;
           case 's':
             outputdepth = CV_8U;
             break;
           case 'i':
             outputdepth = CV_16U;
             break;
           case 'g':
             printf("Forcing grayscale format...\n");
             colformat=CV_LOAD_IMAGE_GRAYSCALE;
             break;
           case 'o':
             printf("Images of different sizes...\n");
             oddsizes=1;
             break;
           default:
             abort ();
           }

    if (argc-optind<7)
    {
        fprintf(stderr, "At least two input images (with offsets) and an output are required.\n");
        return 1;
    }

    nextimgfile = argv[optind];
    
    nextimage = imread( nextimgfile, colformat );
    if ( !nextimage.data )
    {
        printf("Failed to load %s\n",nextimgfile);
        return -1;
    }
    nChannels = (colformat==CV_LOAD_IMAGE_GRAYSCALE) ? 1 : nextimage.channels();

    if (outputdepth==(-1)){
        outputdepth=nextimage.depth();
        outputtype = nextimage.type();
    } else if (outputdepth==CV_8U){
        outputtype=((nChannels==1) ? CV_8UC1 : CV_8UC3);
    } else if (outputdepth==CV_16U){
        outputtype=((nChannels==1) ? CV_16UC1 : CV_16UC3);
    }

    indatamax = ReturnTypeMax( nextimage.depth() );
    outdatamax = ReturnTypeMax( outputdepth );    
    //switch (outputdepth) {
    //    case CV_8U:
    //    case CV_8S:
    //        elementsize = 1L; break;
    //    case CV_16U:
    //    case CV_16S:
    //        elementsize = 2L; break;
    //    case CV_32S:
    //    case CV_32F:
    //        elementsize = 4L; break;
    //    case CV_64F:
    //        elementsize = 8L; break;
    //}  
    //elementsize = 4L; //(outputtype>0) ? outputtype/4L : (2147483648+outputtype)/4L ;
    //imagebytesize = (long int) output_width * (long int) output_height * (long int) nChannels * elementsize;
    //printf("elementsize = %ld, imagebytesize = %ld, INT_MAX=%d\n",elementsize,imagebytesize,INT_MAX); //BJN
    //if (imagebytesize>2147483647L) //can i get rid of this check now?
    //{
    //    printf("Required image size is too large (%ld bytes)...\n",imagebytesize);
    //    return -1;
    //}

    //printf("A...%d\n",INT_MAX); //BJN
    //printf("output_width=%d, output_height=%d, nChannels=%d\n",output_width,output_height,nextimage->nChannels); //BJN  

    overlayimage = Mat(Size(output_width,output_height), ((nChannels==1) ? CV_32FC1 : CV_32FC3) );
    overlayimage.setTo( ((nChannels==1) ? Scalar(0) : Scalar(0,0,0)) );

    tdist = Mat(Size(output_width,output_height), ((nChannels==1) ? CV_32FC1 : CV_32FC3) );
    tdist.setTo( ((nChannels==1) ? Scalar(0) : Scalar(0,0,0)) );

    if (!oddsizes){ //runs faster if you only have to make the distance map once
        cdist = DistMap(nextimage.cols,nextimage.rows,nChannels);
    }

    while (j+optind<argc-1)
    {
        nextimgfile = argv[optind+j];
        offset_x = atof( argv[optind+j+1] );
        offset_y = atof( argv[optind+j+2] );

        nextimage = imread( nextimgfile, colformat );
        if ( !nextimage.data )
        {
            printf("Failed to load %s\n",nextimgfile);
            return -1;
        }
        
        //file conversion
        nextimage.convertTo(next32, CV_32F); 
        if (oddsizes){
            cdist = DistMap(nextimage.cols,nextimage.rows,nChannels);
        }
        multiply(next32,cdist,next32);
        AddROI(tdist, cdist, Point(offset_x,offset_y));
        AddROI(overlayimage, next32, Point(offset_x,offset_y));

        next32.release();

        j=j+3;
    }
    tdist.setTo(1.0,tdist< ((nextimage.rows>nextimage.cols) ? 1.0/(float) nextimage.rows : 1.0/(float) nextimage.cols ));
    divide(overlayimage,tdist,overlayimage);

    overlayimage.convertTo( outimage, outputdepth, outputscale*outdatamax/indatamax );

    overlayimgfile=argv[argc-1];
    imwrite(overlayimgfile,outimage);

}    


