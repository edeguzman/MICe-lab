#include <opencv2/opencv.hpp>

using namespace cv;
using namespace std;

//compile with: 
// g++ cvFilter.cpp -o cvFilter -L/usr/local/lib -I/usr/local/include/opencv -L${exec_prefix}/lib -lopencv_shape -lopencv_stitching -lopencv_objdetect -lopencv_superres -lopencv_videostab -lopencv_calib3d -lopencv_features2d -lopencv_highgui -lopencv_videoio -lopencv_imgcodecs -lopencv_video -lopencv_photo -lopencv_ml -lopencv_imgproc -lopencv_flann -lopencv_viz -lopencv_core


Mat genKernel_d2gauss(double filtwidth,int kernelwidth)
{
    Mat kernel_d2gauss = Mat::zeros(kernelwidth,kernelwidth,CV_64F);
    int j,k;
    double jc,kc,exparg;
    double ksum;
    for (j=0; j<kernel_d2gauss.rows; j++){
        jc = (double) (j-kernelwidth/2);
        for (k=0; k<kernel_d2gauss.cols; k++){
            kc = (double) (k-kernelwidth/2);
            exparg = (pow(jc,2)+pow(kc,2))/(2.0*pow(filtwidth,2));
            kernel_d2gauss.at<double>(Point(k,j)) = (exparg-1.0)*exp(-exparg);
            ksum += kernel_d2gauss.at<double>(Point(k,j)); 
        }
    }
    add(kernel_d2gauss,ksum,kernel_d2gauss);
    kernel_d2gauss /= (-1.0*kernel_d2gauss.at<double>(Point(kernelwidth/2,kernelwidth/2)));
    return kernel_d2gauss;
}

int kernelwidth;
double filtwidth;

//char* outputfile,kernelstr;
string outputfile,kernelstr;

Mat inimage,kernel,filtimage,dximage,dyimage;

int main( int argc, char** argv )
{
    int c;
    char *filename;
     
    while ((c = getopt (argc, argv, "k:w:s:")) != -1)
         switch (c)
         {
           case 'k':
             kernelstr = optarg;
             break;
           case 'w':
             kernelwidth = atoi(optarg);
             break;
           case 's':
             filtwidth = atof(optarg);
             break;
        }

    inimage = imread( argv[optind],CV_LOAD_IMAGE_UNCHANGED);
    if( !inimage.data )
    {
        printf("Failed to load %s\n",argv[optind]);
        return -1;
    }

    outputfile = argv[argc-1];
    
    //setup kernel
    if (kernelstr.compare(0,5,"gauss")==0){ 
        GaussianBlur(inimage,filtimage,Size(kernelwidth,kernelwidth),filtwidth,filtwidth,BORDER_DEFAULT);
    } else if (kernelstr.compare(0,7,"laplace")==0){
        if (filtwidth>0.001){
            GaussianBlur(inimage,inimage,Size(kernelwidth,kernelwidth),filtwidth,filtwidth,BORDER_DEFAULT);
        }
        Laplacian(inimage,filtimage,CV_32F,kernelwidth,1.0,0.0,BORDER_DEFAULT);
        filtimage.convertTo(filtimage,inimage.depth());
    } else if (kernelstr.compare(0,6,"median")==0){
        medianBlur(inimage,filtimage,kernelwidth); //kernelwidth=3||5 supports CV_8U,16U,32F; >5 only supports 8U
    } else if (kernelstr.compare(0,6,"scharr")==0){
        if (filtwidth>0.001){
            GaussianBlur(inimage,inimage,Size(kernelwidth,kernelwidth),filtwidth,filtwidth,BORDER_DEFAULT);
        }
        if (kernelstr.compare(0,7,"scharrx")==0)
            Scharr(inimage,filtimage,CV_32F,1,0,1.0,0.0,BORDER_DEFAULT);
        else if (kernelstr.compare(0,7,"scharry")==0)
            Scharr(inimage,filtimage,CV_32F,0,1,1.0,0.0,BORDER_DEFAULT);
        else {
            Scharr(inimage,dximage,CV_32F,1,0,1.0,0.0,BORDER_DEFAULT);
            Scharr(inimage,dyimage,CV_32F,0,1,1.0,0.0,BORDER_DEFAULT);
            pow(dximage,2.0,dximage);
            pow(dyimage,2.0,dyimage);
            add(dximage,dyimage,filtimage);
            sqrt(filtimage,filtimage); 
        }
        filtimage.convertTo(filtimage,inimage.depth());
    } else if (kernelstr.compare(0,7,"d2gauss")==0){
        Mat kernel_d2gauss = genKernel_d2gauss(filtwidth,kernelwidth);
        filter2D(inimage,filtimage,CV_32F,kernel_d2gauss,Point(-1,-1),0.0,BORDER_DEFAULT);
        filtimage.convertTo(filtimage,inimage.depth()); 
    }
        
    imwrite(outputfile,filtimage);

    return 0;
}
