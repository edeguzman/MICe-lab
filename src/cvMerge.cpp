#include <opencv2/opencv.hpp>
using namespace cv;
using namespace std;

//compile with: 
// g++ cvMerge.cpp -o cvMerge -L/usr/local/lib -I/usr/local/include/opencv -L${exec_prefix}/lib -lopencv_shape -lopencv_stitching -lopencv_objdetect -lopencv_superres -lopencv_videostab -lopencv_calib3d -lopencv_features2d -lopencv_highgui -lopencv_videoio -lopencv_imgcodecs -lopencv_video -lopencv_photo -lopencv_ml -lopencv_imgproc -lopencv_flann -lopencv_viz -lopencv_core

char* outputfile;

Mat red,green,blue;
Mat cred,cgreen,cblue,outimg;
vector<Mat> imgchans; //BGR ordering in channels

int Rin=0,Gin=0,Bin=0;;
int intype;

int main( int argc, char** argv )
{
    int c;
    char *filename;
    opterr = 0;
     
    while ((c = getopt (argc, argv, "r:g:b:")) != -1)
         switch (c)
         {
           case 'r':
             red = imread(optarg,CV_LOAD_IMAGE_UNCHANGED); //CV_LOAD_IMAGE_GRAYSCALE);
             Rin = 1; 
             break;
           case 'g':
             green = imread(optarg,CV_LOAD_IMAGE_UNCHANGED); //CV_LOAD_IMAGE_GRAYSCALE);
             Gin = 1;  
             break;
           case 'b':
             blue = imread(optarg,CV_LOAD_IMAGE_UNCHANGED); //CV_LOAD_IMAGE_GRAYSCALE);
             Bin = 1; 
             break;
         }
    outputfile=argv[argc-1];

    if (Rin==0){
        //printf("No red channel specified...\n");
        if (Gin!=0)
            red = Mat::zeros(green.rows,green.cols,green.type());
        else
            red = Mat::zeros(blue.rows,blue.cols,blue.type());
    }
    if (Gin==0){
        //printf("No green channel specified...\n");
        if (Rin!=0)
            green = Mat::zeros(red.rows,red.cols,red.type());
        else
            green = Mat::zeros(blue.rows,blue.cols,blue.type());
    }
    if (Bin==0){
        //printf("No blue channel specified...\n");
        if (Rin!=0)
            blue = Mat::zeros(red.rows,red.cols,red.type());
        else
            blue = Mat::zeros(green.rows,green.cols,green.type());
    }

    intype = max(max(red.type(),green.type()),blue.type());
 
    if (red.type()==intype)
        cred=red;
    else
        red.convertTo(cred,intype);
    if (green.type()==intype)
        cgreen=green;
    else
        green.convertTo(cgreen,intype);
    if (blue.type()==intype)
        cblue=blue;
    else
        blue.convertTo(cblue,intype);

    imgchans.push_back(blue);
    imgchans.push_back(green);
    imgchans.push_back(red);

    //outimg = Mat::zeros(red.rows,red.cols,CV_16UC3);
    merge(imgchans,outimg);

    imwrite(outputfile,outimg);

}    


