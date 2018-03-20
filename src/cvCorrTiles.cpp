#include <opencv2/opencv.hpp>
#include <glob.h>
using namespace cv;
using namespace std;


//compile with: 
// g++ cvCorrTiles.cpp -o cvCorrTiles -L/usr/local/lib -I/usr/local/include/opencv -L${exec_prefix}/lib -lopencv_shape -lopencv_stitching -lopencv_objdetect -lopencv_superres -lopencv_videostab -lopencv_calib3d -lopencv_features2d -lopencv_highgui -lopencv_videoio -lopencv_imgcodecs -lopencv_video -lopencv_photo -lopencv_ml -lopencv_imgproc -lopencv_flann -lopencv_viz -lopencv_core

void glob_to_strvec(char* globpat,vector<String> &retgloblist)
{
    glob_t glob_result;
    int j;
    glob(globpat,GLOB_TILDE,NULL,&glob_result);
    for (j=0; j<glob_result.gl_pathc; j++){
        retgloblist.push_back(String(glob_result.gl_pathv[j]));
    }
    globfree(&glob_result);
}


char* avgfile;
String postfix_corr_file ("_avgcorr");
String outfile;
Mat cimg,cimg32F,coutimg;
Mat avgimg,avgimg32F,avgimg32FC;
Scalar avgimg_mean;
int j,globinput=0,verbose=0,nChannels,median_kernel_size;
vector<String> imagenamelist; 
vector<String> globvec;
vector<Mat> channels(3);

int main( int argc, char** argv )
{
    int c;
    char *filename;
    opterr = 0;
     
    while ((c = getopt (argc, argv, "gva:m:p:")) != -1)
         switch (c)
         {
           case 'g':
             globinput = 1;
             break;
           case 'v':
             verbose = 1;
             break;
           case 'a':
             avgfile = optarg;
             break;
           case 'm':
             median_kernel_size = atoi(optarg);
             break;
           case 'p':
             postfix_corr_file = String(optarg);
             break;
         }

    //generate input file list
    for (j=optind; j<argc; j++){
        if (globinput){
            glob_to_strvec(argv[j],imagenamelist); //globvec);
        } else {
            imagenamelist.push_back(argv[j]);
        }
    }

    //median filter average image to ensure it is smooth
    avgimg = imread( avgfile, CV_LOAD_IMAGE_UNCHANGED );
    medianBlur(avgimg,avgimg,median_kernel_size);

    //normalize avgimg so it is near unity
    avgimg.convertTo(avgimg32F,CV_32F);
    avgimg_mean=mean(avgimg32F);
    if (avgimg.channels()==1){
        avgimg32F = avgimg32F*(1.0/avgimg_mean.val[0]); 
    } else {
        split(avgimg32F,channels);
        for (j=0; j<avgimg.channels(); j++){
            channels[j]=channels[j]*(1.0/avgimg_mean.val[j]);
        }
        merge(channels,avgimg32F);
    }

    //adapt avgimg32F to same number channels as input
    cimg = imread( imagenamelist[0], CV_LOAD_IMAGE_UNCHANGED );
    nChannels = cimg.channels();
    if (avgimg.channels()<nChannels){
       cvtColor(avgimg32F,avgimg32FC,CV_GRAY2BGR); 
    } else {
       avgimg32FC=avgimg32F;
    }
    threshold(-1.0*avgimg32FC,avgimg32FC,0.1,-1,THRESH_TRUNC);
    avgimg32FC *= -1.0;

    //for each image, multiply tile image by modified avg
    for (j=0; j<imagenamelist.size(); j++){
        if (verbose){ 
            printf("%s \n",imagenamelist[j].c_str()); 
        } else {
            fprintf(stderr,".");
        }
        cimg = imread( imagenamelist[j], CV_LOAD_IMAGE_UNCHANGED );
        cimg.convertTo(cimg32F,((nChannels==1) ? CV_32FC1 : CV_32FC3)); 
        cimg.release();
        divide(cimg32F,avgimg32FC,cimg32F,1.0);
        cimg32F.convertTo(coutimg,cimg.type());
        outfile = imagenamelist[j].substr(0,imagenamelist[j].length()-4) + \
                  String(postfix_corr_file) + \
                  imagenamelist[j].substr(imagenamelist[j].length()-4); 
        imwrite(outfile,coutimg);
        cimg.release();
        cimg32F.release();
        coutimg.release();
    }
    fprintf(stderr,"\n"); 
}    


