
/**
 * Code for thinning a binary image using Zhang-Suen algorithm.
 */

//compile with: gcc -o CGskeletonlength -L/usr/local/lib -I/usr/local/include/opencv CGskeletonlength.c -lopencv_core -lopencv_highgui
//              gcc -o CGskeletonlength CGskeletonlength.c -L/usr/local/lib -I/usr/local/include/opencv -lopencv_core -lopencv_highgui `pkg-config --cflags --libs opencv`


#ifdef _CH_
#pragma package <opencv>
#endif

#define CV_NO_BACKWARD_COMPATIBILITY

#ifndef _EiC
#include "cv.h"
#endif

#include "stdio.h"
#include "unistd.h"

/* 
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
*/

/**
 * Perform one thinning iteration.
 * Normally you wouldn't call this function directly from your code.
 *
 * @param  im    Binary image with range = 0-1
 * @param  iter  0=even, 1=odd
 */
void assignPixelLength(IplImage* im, IplImage* pixlengthimg)
{
    int p2,p3,p4,p5,p6,p7,p8,p9;
    int i,j;
    float val;
    
    for (i = 1; i < im->height-1; i++)
    {
        for (j = 1; j < im->width-1; j++)
        {
            if (im->imageData[im->widthStep * i + j ]==0) continue;
            p2 = (int) im->imageData[im->widthStep * (i-1) + (j) ];   //rem *3 on j terms
            p3 = (int) im->imageData[im->widthStep * (i-1) + (j+1) ]; //im.at<uchar>(i-1, j+1);
            p4 = (int) im->imageData[im->widthStep * (i) + (j+1) ]; //im.at<uchar>(i, j+1);
            p5 = (int) im->imageData[im->widthStep * (i+1) + (j+1) ]; //im.at<uchar>(i+1, j+1);
            p6 = (int) im->imageData[im->widthStep * (i+1) + (j) ]; //im.at<uchar>(i+1, j);
            p7 = (int) im->imageData[im->widthStep * (i+1) + (j-1) ]; //im.at<uchar>(i+1, j-1);
            p8 = (int) im->imageData[im->widthStep * (i) + (j-1) ]; //im.at<uchar>(i, j-1);
            p9 = (int) im->imageData[im->widthStep * (i-1) + (j-1) ]; //im.at<uchar>(i-1, j-1);
            val = 0.5*((float) (p2+p4+p6+p8)) + 0.5*sqrt(2)*((float) (p3+p5+p7+p9));
            pixlengthimg->imageData[pixlengthimg->widthStep * i +j] = (uchar) ( 0.5 + 255.0*val/(2+2*sqrt(2)) ); //0.5 for rounding
        }
    }
}


/**
 * This is an example on how to call the thinning function above.
 */
int main( int argc, char** argv )
{
    int c;
    char *origimgfile;
    char *outputfile;
    opterr = 0;
    IplImage *src=0, *srcgray=0, *bw=0, *pixlengthimg;
    
    while ((c = getopt (argc, argv, "t:m:h")) != -1)
    {
         switch (c)
           {
           case 'h':
             printf("CGskeletonlength [options] <inputimg> <outputimg>\n");
             return -1;
             break;
           default:
             abort ();
           }
    }
    
    if (argc-optind<2)
    {
        fprintf(stderr, "An input and an output image are required arguments.\n");
        return 1;
    }
    origimgfile = argv[optind];
    outputfile = argv[optind+1];

    if( (src = cvLoadImage( origimgfile, 1)) == 0 )
    {
        printf("Failed to load %s\n",origimgfile);
        return -1;
    }
    srcgray = cvCreateImage(cvGetSize(src), IPL_DEPTH_8U, 1);
    cvCvtColor(src, srcgray, CV_RGB2GRAY);

    bw = cvCreateImage(cvGetSize(srcgray), IPL_DEPTH_8U, 1);
    cvThreshold(srcgray, bw, 1, 1, CV_THRESH_BINARY); 

    //output to file
    pixlengthimg = cvCreateImage(cvSize(bw->width,bw->height), IPL_DEPTH_8U, 1); 
    cvSetZero(pixlengthimg);
    assignPixelLength(bw,pixlengthimg);
    cvSaveImage(outputfile,pixlengthimg);
    
    return 0;
}