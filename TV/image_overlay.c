#ifdef _CH_
#pragma package <opencv>
#endif

#define CV_NO_BACKWARD_COMPATIBILITY

#ifndef _EiC
#include "cv.h"
#endif

#include "stdio.h"
#include "unistd.h"

//compile with: gcc -o image_overlay -L/usr/local/lib -I/usr/local/include/opencv image_overlay.c -lopencv_core -lopencv_highgui
//              gcc -o image_overlay image_overlay.c -L/usr/local/lib -I/usr/local/include/opencv -lopencv_core -lopencv_highgui `pkg-config --cflags --libs opencv`

#define CV_LOAD_IMAGE_COLOR       1
#define CV_LOAD_IMAGE_GRAYSCALE   0
#define CV_LOAD_IMAGE_UNCHANGED  -1

#define NONETYPE -1

void AddROI(IplImage* src, IplImage* overlay, CvPoint location)
{
    int x,y,i;
    int ROIwidth,ROIheight;
    ROIwidth = (location.x+overlay->width > src->width) ? src->width-location.x : overlay->width;
    ROIheight = (location.y+overlay->height > src->height) ? src->height-location.y : overlay->height;
    if ((ROIwidth>0)&&(ROIheight>0))
    {
        cvSetImageROI(src,cvRect( location.x, location.y, ROIwidth, ROIheight ));
        cvSetImageROI(overlay,cvRect( 0,0,ROIwidth,ROIheight));
        cvAdd(src,overlay,src,NULL);
        cvResetImageROI(src);
        cvResetImageROI(overlay);
    }
}

IplImage* DistMap(int width, int height, int nChannels)
{
    int x,y,midy,midx;
    int ROIwidth,ROIheight;
    float cdist;
    IplImage* result;
    result = cvCreateImage( cvSize(width, height), IPL_DEPTH_32F, nChannels );
    midx = width/2;
    midy = height/2;
    for (x=0; x<width; x++)
    {
        for (y=0; y<height; y++)
        {
            if (y<midy)
                cdist=(float) y; 
            else
                cdist=(float) (height-y-1); 
            if (x<midx)
                cdist=((float) x > cdist) ? cdist : 1.0+(float) x; 
            else
                cdist=((float) (width-x-1) > cdist) ? cdist : (float) (width-x-1); 
            cvSet2D(result, y, x, cvScalar(cdist,cdist,cdist,cdist));
        }
    }
    return result;
}

int output_width=0,output_height=0,j;
int outputtype=NONETYPE;
int colformat=CV_LOAD_IMAGE_UNCHANGED;
int nChannels=1;
long int elementsize=4L;
long int imagebytesize;

char* nextimgfile;
char* overlayimgfile;

double minVal,maxVal;
double offset_x,offset_y;
double datamax;
double outputscale=1.0;

IplImage *nextimage=0, *overlayimage=0;
IplImage *next32=0, *mov32=0, *out32=0, *outimage=0;
IplImage *cdist = 0, *tdist = 0;
CvPoint offset;
CvPoint minLoc, maxLoc;


int main( int argc, char** argv )
{
    int c;
    char *filename;
    opterr = 0;
     
    while ((c = getopt (argc, argv, "x:y:I:sig")) != -1)
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
             outputtype = IPL_DEPTH_8U;
             break;
           case 'i':
             outputtype = IPL_DEPTH_16U;
             break;
           case 'g':
             printf("Forcing grayscale format...\n");
             colformat=CV_LOAD_IMAGE_GRAYSCALE;
           default:
             abort ();
           }

    if (argc-optind<7)
    {
        fprintf(stderr, "At least two input images (with offsets) and an output are required.\n");
        return 1;
    }

    nextimgfile = argv[optind];
    if( (nextimage = cvLoadImage( nextimgfile, CV_LOAD_IMAGE_UNCHANGED)) == 0 )
    {
        printf("Failed to load %s\n",nextimgfile);
        return -1;
    }
    switch (nextimage->depth)
    {
        case IPL_DEPTH_32F:
        case IPL_DEPTH_64F:
            datamax=1.0;          break;
        case IPL_DEPTH_8U:
            datamax=255.0;        break;
        case IPL_DEPTH_8S:
            datamax=127.0;        break;
        case IPL_DEPTH_16U:
            datamax=65535.0;      break;
        case IPL_DEPTH_16S:
            datamax=32767.0;      break;
        case IPL_DEPTH_32S:
            datamax=2147483647.0; break;
    }
    if (outputtype==NONETYPE) outputtype=nextimage->depth;
    
    nChannels = (colformat==CV_LOAD_IMAGE_GRAYSCALE) ? 1 : nextimage->nChannels;
    elementsize = 4L; //(outputtype>0) ? outputtype/4L : (2147483648+outputtype)/4L ;
    imagebytesize = (long int) output_width * (long int) output_height * (long int) nChannels * elementsize;
    //printf("elementsize = %ld, imagebytesize = %ld, INT_MAX=%d\n",elementsize,imagebytesize,INT_MAX); //BJN
    if (imagebytesize>2147483647L)
    {
        printf("Required image size is too large (%ld bytes)...\n",imagebytesize);
        return -1;
    }

    //printf("A...%d\n",INT_MAX); //BJN
    //printf("output_width=%d, output_height=%d, nChannels=%d\n",output_width,output_height,nextimage->nChannels); //BJN  

    overlayimage = cvCreateImage(cvSize(output_width,output_height), IPL_DEPTH_32F, nextimage->nChannels);
    cvSet(overlayimage, cvScalar(0,0,0,0), NULL);

    tdist = cvCreateImage(cvSize(output_width,output_height), IPL_DEPTH_32F, nChannels);
    cvSet(tdist, cvScalar(0.01,0.01,0.01,0.01), NULL);

    next32 = cvCreateImage(cvGetSize(nextimage), IPL_DEPTH_32F, nChannels);
    cdist = DistMap(next32->width,next32->height,nChannels);

    j=0;
    while (j+optind<argc-1)
    {
        nextimgfile = argv[optind+j];
        offset_x = atof( argv[optind+j+1] );
        offset_y = atof( argv[optind+j+2] );

        if( (nextimage = cvLoadImage( nextimgfile, colformat)) == 0 )
        {
            printf("Failed to load %s\n",nextimgfile);
            return -1;
        }

        //file conversion
        next32 = cvCreateImage(cvGetSize(nextimage), IPL_DEPTH_32F, nChannels);
        cvConvertScale(nextimage, next32, 1.0/datamax, 0.0);

        //cdist = DistMap(next32->width,next32->height,next32->nChannels);-->I moved this outside the loop, assuming all images are the same size 
        cvMul(next32,cdist,next32,1.0);
        AddROI(tdist, cdist, cvPoint(offset_x,offset_y));
        AddROI(overlayimage, next32, cvPoint(offset_x,offset_y));

        cvReleaseImage(&next32);

        j=j+3;
    }

    cvDiv(overlayimage,tdist,overlayimage,1.0);
    out32=overlayimage;
    //cvSetImageCOI(out32, 1);
    //cvMinMaxLoc(out32,&minVal,&maxVal,&minLoc,&maxLoc,NULL);
    //cvSetImageCOI(out32, 0);

    switch (outputtype)
    {
        case IPL_DEPTH_32F:
        case IPL_DEPTH_64F:
            datamax=1.0;          break;
        case IPL_DEPTH_8U:
            datamax=255.0;        break;
        case IPL_DEPTH_8S:
            datamax=127.0;        break;
        case IPL_DEPTH_16U:
            datamax=65535.0;      break;
        case IPL_DEPTH_16S:
            datamax=32767.0;      break;
        case IPL_DEPTH_32S:
            datamax=2147483648.0; break;
    }
    datamax = datamax*outputscale;

    outimage = cvCreateImage(cvGetSize(out32), outputtype, nChannels);

    cvConvertScale(out32, outimage, datamax, 0.0);

    overlayimgfile=argv[argc-1];
    cvSaveImage(overlayimgfile,outimage);

}    


