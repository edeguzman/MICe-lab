#include <opencv2/opencv.hpp>

using namespace cv;
using namespace std;

//compile with: 
// g++ cvRectCrop.cpp -o cvRectCrop -L/usr/local/lib -I/usr/local/include/opencv -L${exec_prefix}/lib -lopencv_shape -lopencv_stitching -lopencv_objdetect -lopencv_superres -lopencv_videostab -lopencv_calib3d -lopencv_features2d -lopencv_highgui -lopencv_videoio -lopencv_imgcodecs -lopencv_video -lopencv_photo -lopencv_ml -lopencv_imgproc -lopencv_flann -lopencv_viz -lopencv_core

int crpleft=0,crptop=0,crpbottom=0,crpright=0;
int xp,yp,xwidth,yheight;

char* outputfile;

Mat inimage,crpimage;

int main( int argc, char** argv )
{
    int c;
    char *filename;
     
    while ((c = getopt (argc, argv, "l:r:t:b:i:")) != -1)
         switch (c)
         {
           case 'l':
             crpleft = atoi(optarg);
             break;
           case 'r':
             crpright = atoi(optarg);
             break;
           case 't':
             crptop = atoi(optarg);
             break;
           case 'b':
             crpbottom = atoi(optarg);
             break;
           case 'i':
             crpleft = crpright = crptop = crpbottom = atoi(optarg);
             break;
         }

    inimage = imread( argv[optind],CV_LOAD_IMAGE_UNCHANGED);
    if( !inimage.data )
    {
        printf("Failed to load %s\n",argv[optind]);
        return -1;
    }

    outputfile = argv[argc-1];

    xp = crpleft; 
    yp = crptop;
    xwidth = inimage.rows-xp-crpright;
    yheight = inimage.cols-yp-crpbottom;

    if ((xwidth<=0)||(yheight<=0)){
        printf("Specified crop produces no image!\n");
        abort();
    }

    crpimage = inimage(Rect(xp,yp,xwidth,yheight));
    imwrite(outputfile,crpimage);

    return 0;
}
