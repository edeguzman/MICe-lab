#include <opencv2/opencv.hpp>
#include <glob.h>
using namespace cv;
using namespace std;


//compile with: 
// g++ cvAvgImages.cpp -o cvAvgImages -L/usr/local/lib -I/usr/local/include/opencv -L${exec_prefix}/lib -lopencv_shape -lopencv_stitching -lopencv_objdetect -lopencv_superres -lopencv_videostab -lopencv_calib3d -lopencv_features2d -lopencv_highgui -lopencv_videoio -lopencv_imgcodecs -lopencv_video -lopencv_photo -lopencv_ml -lopencv_imgproc -lopencv_flann -lopencv_viz -lopencv_core

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


char* outputfile;
float Iscale;
Mat cimg,cimg32F;
Mat avgimg;
int N,j,globinput=0,verbose=0,maximgs=-1,nChannels,outtype;
vector<String> imagenamelist; 
vector<String> globvec;

int main( int argc, char** argv )
{
    int c;
    char *filename;
    opterr = 0;
     
    while ((c = getopt (argc, argv, "gvs:m:")) != -1)
         switch (c)
         {
           case 'g':
             globinput = 1;
             break;
           case 'v':
             verbose = 1;
             break;
           case 's':
             Iscale = atof(optarg);
             break;
           case 'm':
             maximgs = atoi(optarg);
             break;
         }
    outputfile=argv[argc-1];

    //generate input file list
    for (j=optind; j<argc-1; j++){
        printf("%s\n",argv[j]);
        if (globinput){
            glob_to_strvec(argv[j],imagenamelist); //globvec);
        } else {
            imagenamelist.push_back(argv[j]);
        }
    }

    //make average accumulator
    cimg = imread( imagenamelist[0], CV_LOAD_IMAGE_UNCHANGED );
    nChannels = cimg.channels();
    outtype = cimg.type();
    avgimg = Mat(Size(cimg.cols,cimg.rows), ((nChannels==1) ? CV_32FC1 : CV_32FC3) );
    avgimg.setTo( ((nChannels==1) ? Scalar(0) : Scalar(0,0,0)) );

    //add image one at a time
    N=imagenamelist.size();
    if (maximgs>0) N = (N>maximgs) ? maximgs : N;
    printf("N = %d\n",N);
    for (j=0; j<N; j++){
        fprintf(stderr,".");
        if (verbose){ printf("%s \n",imagenamelist[j].c_str()); }
        cimg = imread( imagenamelist[j], CV_LOAD_IMAGE_UNCHANGED );
        cimg.convertTo(cimg32F,((nChannels==1) ? CV_32FC1 : CV_32FC3)); //scaling here causes discretization problems, do it at accumulate step
        cimg.release();
        accumulate(cimg32F*Iscale/((float) N),avgimg);
        cimg32F.release();
    }
    fprintf(stderr,"Done Averaging\n");
    avgimg.convertTo(avgimg,outtype);
    
    //output image
    imwrite(outputfile,avgimg);

}    


