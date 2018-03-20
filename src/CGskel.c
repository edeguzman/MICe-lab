
/**
 * Code for thinning a binary image using Zhang-Suen algorithm.
 */

//compile with: gcc -o CGskel -L/usr/local/lib -I/usr/local/include/opencv CGskel.c -lopencv_core -lopencv_highgui
//              gcc -o CGskel CGskel.c -L/usr/local/lib -I/usr/local/include/opencv -lopencv_core -lopencv_highgui `pkg-config --cflags --libs opencv`


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
void thinningIteration(IplImage* im, int iter)
{
    IplImage* marker = cvCreateImage(cvSize(im->width,im->height), IPL_DEPTH_8U, 1); 
    cvSetZero(marker);
    uchar p2,p3,p4,p5,p6,p7,p8,p9;
    int i,j;
    
    for (i = 1; i < im->height-1; i++)
    {
        for (j = 1; j < im->width-1; j++)
        {
            p2 = im->imageData[im->widthStep * (i-1) + (j) ];   //rem *3 on j terms
            p3 = im->imageData[im->widthStep * (i-1) + (j+1) ]; //im.at<uchar>(i-1, j+1);
            p4 = im->imageData[im->widthStep * (i) + (j+1) ]; //im.at<uchar>(i, j+1);
            p5 = im->imageData[im->widthStep * (i+1) + (j+1) ]; //im.at<uchar>(i+1, j+1);
            p6 = im->imageData[im->widthStep * (i+1) + (j) ]; //im.at<uchar>(i+1, j);
            p7 = im->imageData[im->widthStep * (i+1) + (j-1) ]; //im.at<uchar>(i+1, j-1);
            p8 = im->imageData[im->widthStep * (i) + (j-1) ]; //im.at<uchar>(i, j-1);
            p9 = im->imageData[im->widthStep * (i-1) + (j-1) ]; //im.at<uchar>(i-1, j-1);

            int A  = (p2 == 0 && p3 == 1) + (p3 == 0 && p4 == 1) + 
                     (p4 == 0 && p5 == 1) + (p5 == 0 && p6 == 1) + 
                     (p6 == 0 && p7 == 1) + (p7 == 0 && p8 == 1) +
                     (p8 == 0 && p9 == 1) + (p9 == 0 && p2 == 1);
            int B  = p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9;
            int m1 = (iter == 0) ? (p2 * p4 * p6) : (p2 * p4 * p8);
            int m2 = (iter == 0) ? (p4 * p6 * p8) : (p2 * p6 * p8);

            if (A == 1 && (B >= 2 && B <= 6) && m1 == 0 && m2 == 0)
                marker->imageData[marker->widthStep * i + j] = 1;
        }
    }
    cvNot(marker,marker);
    cvAnd(im,marker,im,NULL); //im &= ~marker;
}

/**
 * Function for thinning the given binary image
 *
 * @param  im  Binary image with range = 0-255
 */
void thinning(IplImage* im)
{
    double minVal,maxVal;
    cvMinMaxLoc(im,&minVal,&maxVal,NULL,NULL,NULL);
    cvConvertScale(im,im,1.0/maxVal,0.0);

    IplImage* prev = cvCreateImage(cvSize(im->width,im->height), IPL_DEPTH_8U, 1); 
    cvSetZero(prev);
    IplImage* diff = cvCreateImage(cvSize(im->width,im->height), IPL_DEPTH_8U, 1);

    do {
        thinningIteration(im, 0);
        thinningIteration(im, 1);
        cvAbsDiff(im, prev, diff);
        cvCopy(im,prev,NULL);
    } 
    while (cvCountNonZero(diff) > 0);

    cvConvertScale(im,im,maxVal,0.0);
}

/**
 * This is an example on how to call the thinning function above.
 */
int main( int argc, char** argv )
{
    int c;
    int Ithreshold=15,nmorphDE=2;
    char *origimgfile;
    char *skelimgfile;
    opterr = 0;
    IplImage *src=0, *srcgray=0, *bw=0;
    
    while ((c = getopt (argc, argv, "t:m:h")) != -1)
    {
         switch (c)
           {
           case 't':
             Ithreshold = atoi(optarg);
             break;
           case 'm':
             nmorphDE = atoi(optarg);
             break;
           case 'h':
             printf("CGskel [options] <inputimg> <outputimg>\n");
             printf("    -t <int> threshold for binarization (default: %d, assuming uchar 255 max)\n",Ithreshold);
             printf("    -m <int> number of successive dilations followed by successive erosions prior to skeletonization (default %d)\n",nmorphDE);
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
    skelimgfile = argv[optind+1];

    if( (src = cvLoadImage( origimgfile, 1)) == 0 )
    {
        printf("Failed to load %s\n",origimgfile);
        return -1;
    }
    srcgray = cvCreateImage(cvGetSize(src), IPL_DEPTH_8U, 1);
    cvCvtColor(src, srcgray, CV_RGB2GRAY);

    bw = cvCreateImage(cvGetSize(srcgray), IPL_DEPTH_8U, 1);
    cvThreshold(srcgray, bw, Ithreshold, 255, CV_THRESH_BINARY); 

    //dilate N, erode N
    if (nmorphDE>0){
       cvDilate(bw,bw,NULL,nmorphDE);
       cvErode(bw,bw,NULL,nmorphDE);
    }

    //open operation to eliminate isolated pixels
    cvErode(bw,bw,NULL,1);
    cvDilate(bw,bw,NULL,1); 

    //skeletonize
    thinning(bw);

    //output to file
    cvSaveImage(skelimgfile,bw);
    
    return 0;
}