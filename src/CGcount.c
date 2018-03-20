
/**
 * Code for counting unique pixel values in a specified radius around a pixel
 */

//compile with: gcc -o CGcount CGcount.c -L/usr/local/lib -I/usr/local/include/opencv -lopencv_core -lopencv_highgui
//              gcc -o CGcount CGcount.c -L/usr/local/lib -I/usr/local/include/opencv -lopencv_core -lopencv_highgui `pkg-config --cflags --libs opencv`


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


void numunique(IplImage* in, IplImage* out, int radius)
{
    int sqr_rad = radius*radius;
    if (in->nChannels > 1) // || in.type() != CV_32F) 
    {
        printf("numunique only works with 1-channel data\n");
        printf("    image provided has %d channels. \n",in->nChannels);
        return;
    }

    int val_list[radius*radius];
    int y,x,i,j,k,q,clen=0,cval;
    char vflag=0;
    int nbytes = in->depth/8;
    
    for (y = 0; y < in->height; y++)
    {
        for (x = 0; x < in->width; x++)
        {
            clen=0;
            for (i = -radius; i<=radius; i++){
                for (j = -radius; j<=radius; j++){
                    if ( ( (i*i+j*j)>sqr_rad )||
                         ( (y+i)<0 )||
                         ( (y+i)>=in->height )||
                         ( (x+j)<0 )||
                         ( (x+j)>=in->width ) ) continue;
                    cval=0;
                    for (q=0; q<nbytes; q++){
                        cval += (in->imageData[in->widthStep * (y+i) + (nbytes * (x+j)) + q ])<<(8*(nbytes-1-q));
                    }
                    //check if value is new
                    vflag=0;
                    for (k = 0; k<clen; k++){
                        if (cval==val_list[k]){
			  vflag=1;
			  break;
			}
                    }
                    if (vflag==0){
                        val_list[clen]=cval;
                        clen+=1;
                    }
		}
 	    }
	    //length of array gives number unique values
            out->imageData[out->widthStep * (y) + (x) ] = (char) (clen);
        }
    }
}




int main( int argc, char** argv )
{
    int c;
    int radius=21,nmorphDE=2;
    char *origimgfile;
    char *outputfile;
    opterr = 0;
    IplImage *src=0, *outputimg=0;
    
    while ((c = getopt (argc, argv, "r:h")) != -1)
    {
         switch (c)
           {
           case 'r':
             radius = atoi(optarg);
             break;
           case 'h':
             printf("CGcount [options] <inputimg> <outputimg>\n");
             printf("    -r <int> radius of counting ROI (default: %d)\n",radius);
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

    if( (src = cvLoadImage( origimgfile, -1)) == 0 )
    {
        printf("Failed to load %s\n",origimgfile);
        return -1;
    }

    outputimg = cvCreateImage(cvSize(src->width,src->height), IPL_DEPTH_8U, 1);  //need to change depth BJN
    cvSetZero(outputimg);

    //find unique vaiues within radius
    numunique(src,outputimg,radius);

    //output to file
    cvSaveImage(outputfile,outputimg);
    
    return 0;
}