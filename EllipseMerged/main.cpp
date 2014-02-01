 ////////////////////////////////////////////////////////////////////////////////
//                          PROJET				                              //
////////////////////////////////////////////////////////////////////////////////

#include "CImg.h"
#include <iostream>
#include <math.h>
#include <vector>
#include <algorithm>


#define DETECT 20

using namespace cimg_library;
using namespace std;
#define DEGREES_TO_RADIANS(degrees) (degrees*M_PI/180)
#define RADIANS_TO_DEGREES(radians) (radians*180/M_PI)

const int HoughHistoThreshold = 1;
const int HoughKNThreshold = 1;
const int HoughCenterThreshold = 1000;

typedef struct center_point_type {
    int a0, b0;
    center_point_type() = default;
    center_point_type(const int& x, const int& y) : a0(x), b0(y){}
    
} CENTER;

typedef struct angles_type {
    float atanK, atanN;
    angles_type() = default;
    angles_type(const int& angleK, const int& angleN) :
        atanK(angleK), atanN(angleN){}
    
} ANGLES;

typedef struct ellipse_type{
    ellipse_type() = default;
    float a0;
    float b0;
    float atanK;
    float atanN;
    float ax;
    float ay;
    float bx;
    float by;
    float a;
    float b;
    ellipse_type(CENTER c, ANGLES ang, int radioAx) : a0(c.a0), b0(c.b0), atanK(ang.atanK),
        atanN(ang.atanN), ax(radioAx){
        float atanKrad = DEGREES_TO_RADIANS(atanK);
        float atanNrad = DEGREES_TO_RADIANS(atanN);
        float K = tan(atanKrad);
        float N = tan(atanNrad);
        ay = K*ax;
        by = N*ax;
        bx = -ay*by /ax;
        a = sqrt(ax*ax + ay*ay);
        b = sqrt(bx*bx+by*by);
    }
    
} ELLIPSE;


/*******************************************************************************
 
 MaxDetection: Compute local maximum
 ImgIn   : Input Image
 ImgOut  : Maximum map of ImgIn
 
 *******************************************************************************/
void MaxDetection(CImg<> ImgIn, CImg<> &ImgOut)
{
    ImgOut.fill(0);
    CImgList<> ImgInGrad = ImgIn.get_gradient("xy",4);
    
    cimg_forXY(ImgOut,x,y)
    {
        float norm = sqrt(ImgInGrad[0](x,y)*ImgInGrad[0](x,y)+ImgInGrad[1](x,y)*ImgInGrad[1](x,y));
        
        float dx = ImgInGrad[0](x,y)/norm;
        float dy = ImgInGrad[1](x,y)/norm;
        
        if(ImgIn(x,y)>ImgIn.linear_atXY(x+dx,y+dy) && ImgIn(x,y)>ImgIn.linear_atXY(x-dx,y-dy))
                ImgOut(x,y) = ImgIn(x,y);
    }
}

/***************************************************************
  IsLocalMax: Tells if the point (a0,b0) is a local max of an accumulator
 Acc:   Accumulator
 a0 : x coordinate
 b0: y coordinate
 
 *************************************************************/
bool IsLocalMax(CImg<> &Acc, int a0, int b0){
     long max = Acc(a0,b0);
    for(int dx = -4; dx < 4; ++dx)
        for(int dy = -4; dy < 4; ++dy)
            if(a0 + dx < Acc.width() && a0+dx > 0 && b0+dy > 0 && b0+dy < Acc.height()){
                if(Acc(a0+dx, b0+dy) > max)
                    return false;
                
            }
    return true;
}

/***************************************************************
 MaxAccCenterModifyAcc: Computes the maximum of the center accumulator, erasing the point afterwards. 
Returns true if a maximum can be computed, false otherwise.
 Acc:   Accumulator
 a0 : returned x coordinate
 b0: returned y coordinate
 
 *************************************************************/
bool MaxAccCenterModifyAcc(CImg<> &Acc, int& a0, int& b0){
    a0 = -1;
    b0 = -1;
    long max = 0;
    bool isThereACenter = false;
    cimg_forXY(Acc,x,y)
        if(Acc(x,y) > 0 && Acc(x, y) > max){
            if(IsLocalMax(Acc, x , y)){
                max = Acc(x,y);
                a0 = x;
                b0 = y;
            isThereACenter = true;
            }
            
        }
    if(isThereACenter)
        for(int dx = -4; dx < 4; ++dx)
            for(int dy = -4; dy < 4; ++dy)
                if(a0+dx > 0 && a0 + dx < Acc.width() && b0+dy > 0 && b0+dy < Acc.height())
                    Acc(a0+dx, b0+dy) = -1;
                
            
        
    
    
    
    return isThereACenter;
}

/***************************************************************
 proportion_ellipse: How much of an ellipse exists in an image
 ImgIn:   Image
 ellipse : Ellipse
 min_dist: Maximum error allowed
 *************************************************************/
float proportion_ellipse(CImg<>& ImgIn, ELLIPSE ellipse,float min_dist ){
    float count= 0.0;
    float perimeter;
    float a = ellipse.a;
    float b = ellipse.b;
    int x1, y1, x2, y2;
    cimg_forXY(ImgIn, x, y)
        if(ImgIn(x,y) > 0){
            float theta = DEGREES_TO_RADIANS(ellipse.atanK);
            x1 = x - ellipse.a0;
            y1 = y - ellipse.b0;
            x2 = x1*cos(theta) + y1*sin(theta);
            y2 = -x1*sin(theta) + y1*cos(theta);
            
            x2 /= ellipse.a;
            y2 /= ellipse.b;
            if (abs((1.0 - (x2*x2 + y2*y2))) < min_dist)
            {
                count++;
            }
        }
    perimeter = M_PI * ( 3*(a+b) - sqrt((a+3*b)*(3*a+b)) );
    return count/perimeter;
    
    
}

/***************************************************************
 IsLocalMaxOfHisto: Tells if a cell is a local maximum of an histogram
 Acc:   Accumulator
 a0 : x coordinate
 b0: y coordinate
 
 *************************************************************/
bool IsLocalMaxOfHisto(CImg<long>& Histo, int index){
    long max = Histo(index);
    for(int dx = -4; dx < 4; ++dx)
        if(index + dx > 0 && index + dx < Histo.width())
            if(Histo(index + dx) > max)
                return false;
    return true;
            
}

/***************************************************************
 MaxHistoModifyAcc: Compute the cell where the maximum of an histogram is located.
 Histo:   Histogram
 index : Position of the maximum of the histogram
 
 *************************************************************/
bool MaxHistoModifyAcc(CImg<long>& Histo, int& index){
    long max = -1;
    index = -1;
    bool isThereAx = false;
    cimg_forX(Histo, x)
    if(Histo(x) > 0 && Histo(x) > max && Histo(x) < Histo.width() && IsLocalMaxOfHisto(Histo, x)){
        max = Histo(x);
        index  = x;
        isThereAx = true;
    }
    if(isThereAx)
        for(int dx = -4; dx < 4; ++dx)
            if(index + dx > 0 && index + dx < Histo.width())
                Histo(index+dx) = -1;
    return isThereAx;
}

/***************************************************************
 maxAccKNModifyAcc: Computes the position of the maximum value of a KN accumulator
 Acc:   Accumulator
 angles: Position (atanK, atanN) of the maximum of the KN accumulator
 *************************************************************/
bool maxAccKNModifyAcc(CImg<>& Acc, ANGLES& angles){
    int atanK = -1;
    int atanN = -1;
        bool areThereAngles = false;
        long  max = -1;
        for(int x = 0; x < Acc.width(); ++x)
            for(int y = 0; y < Acc.height(); ++y){
                if(Acc(x,y) > 0 && Acc(x,y) > max && IsLocalMax(Acc, x, y)){
                    max = Acc(x,y);
                    atanK = x;
                    atanN = y;
                    
                    areThereAngles = true;
                }
            }
    if(areThereAngles){
        
            for(int dx = -4; dx < 4; ++dx)
                for(int dy = -4; dy < 4; ++dy)
                    if(atanK+dx > 0 && atanK + dx < Acc.width() && atanN+dy > 0 && atanN+dy < Acc.height())
                        Acc(atanK+dx, atanN+dy) = -1;
        
        angles.atanK = atanK;
        angles.atanN = atanN;
}
        
  
    
    return areThereAngles;
    
}

/***************************************************************
 existsEllipse: Tells if a determined ellipse exists in an image.
 ImgIn:   Image
 ellipse: Ellipse to be found
 
 *************************************************************/
bool existsEllipse(CImg<> &ImgIn, const ELLIPSE& ellipse){
    float min_ellipse_proportion = 1.2;
    float a_over_b;
    float a = ellipse.a;
    float b = ellipse.b;
    
    unsigned char purple[] = { 233,100,125 };
    if(a > 0 && b > 0){
        ///////
        float proportion = proportion_ellipse(ImgIn, ellipse, 0.01);
        if(a != 0 && b != 0){
            if(a > b){
                a_over_b = ((float)a) / (float)b;
            } else {
                a_over_b = ((float)b) / ((float)a);
            }
        } else {
            a_over_b = 0;
        }
        //if ( (proportion > min_ellipse_proportion) && (a >= min_ellipse/2.0 && b >= min_ellipse/2.0 ) && a_over_b >= 1.50 )
        if ( (proportion > min_ellipse_proportion)  && a_over_b >= 1.30 )
            return true;

}
    return false;
}

/***************************************************************
 Module: Computes the module of an image
 image:   Image
 *************************************************************/
CImg<> Module (CImg<> image)
{
	//CImgList<float> list = image.get_gradient("xy", 0);
    CImgList<float> G = image.get_gradient("xy", 3);
    CImg<> module = (G[0].sqr()+G[1].sqr()).get_sqrt();
    CImg<> ToReturn(image.width(), image.height());
    cimg_forXY(image, x, y){
        if(module(x,y) > DETECT)
            ToReturn(x,y) = 255;
        else
            ToReturn(x,y) = 0;
    }
	return ToReturn;
    
    
}

/***************************************************************
 AccumulateCenters: Computes the accumulator of centers of ellipses
 ImgIn:   Image
 Acc1: Accumulator of centers
 
 *************************************************************/
bool AccumulateCenters(CImg<> ImgIn, CImg<> &Acc1)
{
    bool isACenterFound = false;
	Acc1.fill(0);
	CImg<> module = Module(ImgIn);
	CImgList<> grad = ImgIn.get_gradient("xy", 3);
    
	cimg_forXY(module, x1, y1)
	{
        if(module(x1,y1) > 0)
        cimg_forXY(module, x2, y2)
        {
            if(module(x2,y2) > 0 && sqrt((y2-y1)*(y2-y1) + (x2-x1)*(x2-x1)) < 25)
                continue;
            
            // Initialization
            float xv = -grad.at(1)(x1, y1);
            float yv = grad.at(0)(x1, y1);
            float xw = -grad.at(1)(x2, y2);
            float yw = grad.at(0)(x2, y2);
            
            float m1 = -xv/yv;
            float m2 = -xw/yw;
            
            if ( module(x1, y1) > DETECT && module(x2, y2) > DETECT && abs(m1-m2) > DEGREES_TO_RADIANS(10)) {
                int xm = (x1+x2)/2, ym = (y1+y2)/2;
                float xt = ((xv*xw)/(yv*xw - yw*xv)) * (y2-y1 + x1*(yv/xv) - x2*(yw/xw));
                float yt = y1 + (xt - x1) * (yv/xv);
                
                if (xt < 0 || xt >= module.width() || yt < 0 || yt >= module.height())
                    continue;
                int b0 = 0;
                
                // Variation de a0
                for (int a0 = 0; a0 < module.width(); ++a0)
                {
                    b0 = (a0 - xm) * ((yt-ym)/(xt-xm)) + ym;
                    
                    
                    if (b0 >= module.height() || b0 < 0)
                        continue;
                    
                    Acc1(a0, b0) += 1;
                    if(!isACenterFound)
                        isACenterFound = true;
                    
                }
                
            }
            
        }
        
    }
    return isACenterFound;
    
}

/***************************************************************
 ExtractCandidateEllipses: Computes all the possible ellipses that could be found in an image
 ImgIn:   Image
 candidateEllipses : Possible ellipses in the image
 centers:  Accumulator of centers to be used
 minA: Minimum length of the small semi-axe
 maxA: Maximum length of the small semi-axe
 minB: Minimum length of the big semi-axe
 maxB: Maximum length of the big semi-axe
 numEllipses: Maximum number of ellipses to be found
 
 *************************************************************/
void ExtractCandidateEllipses(CImg<> &ImgIn, vector<ELLIPSE>& candidateEllipses, const vector<CENTER>& centers, int minA, int maxA, int minB, int maxB, int numEllipses){
    
    int
    xp, yp,atanN, a, b, x0,y0;
    
    float
    slopeFirstDerivativeP, slopeSecondDerivativeP, angleDerivativeP1, angleDerivativeP2, angleFirstDerivativeP, angleSecondDerivativeP, M1, M2, N, K,phi1, phi2, ro, X, Y,ax, ay, bx, by;
    bool isThereACandidate = false;
    CImg<long> Histo(ImgIn.width()/2);
    float margin = 0.0001;
    CImg<> module(ImgIn.width(), ImgIn.height());
    CImg<> Acc2(180, 180);
    module = Module(ImgIn);
    CImgList<> grad = ImgIn.get_gradient("xy", 3);
    for(int i = 0; i < numEllipses; ++i){
        CENTER center = centers.at(i);
            //For each ellipse candidate
        std::cout << "Center: (" << center.a0 << ", " << center.b0 << ")" << std::endl;
        int a0 = center.a0;
        int b0 = center.b0;
            cimg_forXY(module, x1, y1)
                if(module(x1, y1) > 0) //There is an edge point
                    cimg_forXY(module, x2, y2)
                    if(module(x2, y2) > 0 && sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)) > 25 && abs(x2 - x1) > 10 ){
                        // Initialization
                        float xv = -grad.at(1)(x1, y1);
                        float yv = grad.at(0)(x1, y1);
                        float xw = -grad.at(1)(x2, y2);
                        float yw = grad.at(0)(x2, y2);
                        
                        float m1 = -xv/yv;
                        float m2 = -xw/yw;
                        angleDerivativeP1 = atan2(yv, -xv);
                        angleDerivativeP2 = atan2(yw, -xw);
                       
                        
                        if(abs(angleDerivativeP1-DEGREES_TO_RADIANS(90)) < margin || abs(angleDerivativeP1 - DEGREES_TO_RADIANS(270)) < margin || abs(angleDerivativeP2-DEGREES_TO_RADIANS(90)) < margin || abs(angleDerivativeP2 - DEGREES_TO_RADIANS(270)) < margin) continue;
            
                        
                        if(abs(m1-m2) > DEGREES_TO_RADIANS(10)){
                            X = x2 - x1;
                            Y = y2 - y1;
                            
                            
                            if(X==0)
                                continue;
                            slopeFirstDerivativeP = Y/X;
                            float divisor = (X*M2-2*Y);
                            if(divisor <= margin) continue;
                            slopeSecondDerivativeP = (2*M1*X - Y*M2)/divisor;
                            M1 = m1*m2;
                            M2 = m1+m2;
                            int xm = (x1+x2)/2, ym = (y1+y2)/2;
                            
                            float xt = ((xv*xw)/(yv*xw - yw*xv)) * (y2-y1 + x1*(yv/xv) - x2*(yw/xw));
                            float yt = y1 + (xt - x1) * (yv/xv);
                            //Line MT: y - yt = ((ym-yt)/(xm-xt))*(x - xt)
                            //Equation 4: y = slopeSecondDerivativeP*(x - a0) + b0
                            //Find xp -> Find yp
                            float mt = (ym-yt)/(xm-xt);
                            xp = floor(((-mt)*(xt) - b0 + yt + slopeSecondDerivativeP*a0)/(slopeSecondDerivativeP - mt)+0.5f);
                            yp = floor(slopeSecondDerivativeP*(xp - a0) + b0 +0.5f);
                            if(xp >0 && xp < ImgIn.width() && yp > 0 && yp < ImgIn.height()){
                                    phi1 = atan2(Y,X);
                                   
                                    phi2 = atan((2*M1*X - Y*M2)/(X*M2-2*Y));
                                  
                                  
                                
                                    for(float i = 0; i < 180; i+=1){
                                        ro = DEGREES_TO_RADIANS(i);
                                    
                                    if(abs(ro-DEGREES_TO_RADIANS(90)) <= margin || abs(ro-DEGREES_TO_RADIANS(270)) <= margin) continue;
                                        if(i == 0 || i == 90 ) continue;
                                    if(abs((phi1-ro)-DEGREES_TO_RADIANS(90)) <= margin || abs((phi1-ro)-DEGREES_TO_RADIANS(270)) <= margin || abs((phi2-ro) - DEGREES_TO_RADIANS(90)) <= margin || abs((phi2-ro)-DEGREES_TO_RADIANS(270)) <= margin) continue;
                                        N = sqrt(abs(tan(phi1 - ro)*tan(phi2 - ro)));
                                        K = tan(ro);
                                        x0 = floor((xp - a0)/sqrt(K*K+1) + ((yp - b0)*K)/sqrt(K*K+1) + 0.5f);
                                        y0 = floor(((xp - a0)*K)/sqrt(K*K+1) + (yp - b0)/sqrt(K*K+1) +0.5f);
                                        if(sqrt((x0+a0-xp)*(x0+a0-xp)+(y0+b0-yp)*(y0+b0-yp))<= 5){
                                            ax = sqrt(abs((y0*y0 + x0*x0*N*N)/(N*N*(1 + K*K))));
                                            if(ax <= margin)
                                                continue;
                                                    ay = K*ax;
                                                    by = N*ax;
                                                    bx = -ay*by/ax;
                                                    a = sqrt(ax*ax + ay*ay);
                                                    b = sqrt(bx*bx + by*by);
                                                    if(abs(ax + a0) < module.width() && abs(ay + b0) < module.height() && a >= minA && a <= maxA && b >= minB && b <= maxB){
                                                        atanN = RADIANS_TO_DEGREES(atan(N));
                                                        if(atanN < 0)
                                                            atanN +=360;
                                                        
                                                        
                                                        Acc2(floor(i + 0.5f), floor(atanN+0.5f)) += 1;
                                                        std::cout << "Added to Acc2:" << floor(i + 0.5f) << ", " << floor(atanN+0.5f) << std::endl;
                                                        Histo(floor(ax+0.5f)) += 1;
                                                        std::cout << "Added ax: " << floor(ax+0.5f) << std::endl;
                                                        isThereACandidate = true;
                                                        
                                                        
                                                        
                                                        
                                                    }
                                            }
                                    }
                                    
                                    }
                                    }
                                }
    }
        /* Threshold accumulators */
        if(isThereACandidate){
            CImg<> Acc2Thresholded(Acc2.width(), Acc2.height());
            Acc2Thresholded.fill(0);
            cimg_forXY(Acc2, x, y)
                if(Acc2(x,y) < HoughKNThreshold)
                    Acc2(x,y) = 0;
            cimg_forX(Histo, x)
                if(Histo(x) < HoughHistoThreshold)
                    Histo(x) = 0;
                
            
            MaxDetection(Acc2, Acc2Thresholded);
            CImg<> HistoThresholded(Histo);
            HistoThresholded.fill(0);
            MaxDetection(Histo, HistoThresholded);
            vector<ANGLES> list_of_angles;
            
            CImg<> copyAcc2(Acc2Thresholded);
            CImg<long> copyH(HistoThresholded);
            vector<int> list_of_ax;
            for(int i = 0; i < numEllipses; ++i){
                ANGLES angles;
               
                if(maxAccKNModifyAcc(copyAcc2, angles)){
                    list_of_angles.push_back(angles);
                }
                else
                    break;
            }
            for(int i = 0; i < numEllipses; ++i){
                int ax;
                if(MaxHistoModifyAcc(copyH, ax)){
                    list_of_ax.push_back(ax);
                }
                else
                    break;
            }
            
            
            /*
            std::for_each(list_of_angles.begin(), list_of_angles.end(), [&](const ANGLES& theAngles){
                
                    int ax;
                    if(MaxHistoModifyAcc(copyH,ax)){
                        ELLIPSE ellipse(center, theAngles, ax);
                        candidateEllipses.push_back(ellipse);
                    }
            });
             */
            
            //Version 1: With ordering
            /*
            for(int i = 0; i < numEllipses; ++i){
                if(i < centers.size() && i < list_of_angles.size() && i < list_of_ax.size()){
                    ELLIPSE ellipse(centers.at(i), list_of_angles.at(i), list_of_ax.at(i));
                    candidateEllipses.push_back(ellipse);
                }
                else
                    break;
            }
                
             */
            //Version 2:
            std::for_each(centers.begin(), centers.end(), [&](CENTER center){
                std::for_each(list_of_angles.begin(), list_of_angles.end(), [&](ANGLES theAngles){
                    std::for_each(list_of_ax.begin(), list_of_ax.end(), [&](int theIndex){
                        ELLIPSE ellipse(center, theAngles, theIndex);
                        candidateEllipses.push_back(ellipse);
                    });
                });
                
            });
            isThereACandidate = false;
        }
    
}

/*******************************************************************************
 
 Main
 
 *******************************************************************************/
int main(int argc,char **argv)
{
    cimg_usage("Retrieve command line arguments");
    const char* filename = cimg_option("-i","/Users/rubcuevas/Desktop/Algorithmie de l'image/EllipseDetection/EllipseMerged/ellipsefilled.bmp","Input image file");
    
    /* Variables declaration */
    const int minA = 13;
    const int maxA = 17;
    const int minB = 26;
    const int maxB = 33;
    const int numMaxEllipses = 1;
    const unsigned char color[] = {255,255,0};
    CImg<> img(filename);
    CImg<> module = Module(img);
    CImgDisplay moduleSpatial(module, "Module");
    CImgDisplay dispSpatial(img,"Input Image");
    
    /*********************************************************/
    /* Begin accumulation of centers */
    /**********************************************************/
    
    /* Declare and initialize accumulator of centers */
    CImg<> Acc1(img.width(), img.height());
    Acc1.fill(0);
    
    
    /* Accumulate */
    AccumulateCenters(img, Acc1);
    cimg_forXY(Acc1, x, y){
        if(Acc1(x,y) < HoughCenterThreshold)
            Acc1(x,y) = 0;
    }
    /* Threshold centers accumulator */
    CImg<> Acc1Thresholded(Acc1.width(), Acc1.height());
    Acc1Thresholded.fill(0);
    MaxDetection(Acc1, Acc1Thresholded);
    CImg<> Acc1Maxima(Acc1);
    Acc1Maxima.fill(0);

    const CImg<> Acc1Copy(Acc1Thresholded);
    
    /* Display centers accumulator */
    CImgDisplay acc1Spatial(Acc1Copy, "Acc1 Thresholded");
    
    
    /*
      Save centers in a vector of centers.
     */
    vector<CENTER> centers;
    for(int i = 0; i < numMaxEllipses; ++i){
        int a0, b0;
        if(MaxAccCenterModifyAcc(Acc1Thresholded, a0, b0)){
            CENTER c(a0,b0);
            centers.push_back(c);
        }
        
    }
    
    /*********************************************************/
    /* Accumulation of centers finished */
    /**********************************************************/
    vector<ELLIPSE> candidateEllipses;
    ExtractCandidateEllipses(img, candidateEllipses, centers, minA, maxA, minB,maxB,numMaxEllipses);
    int drawn = 0;
    std::for_each(candidateEllipses.begin(), candidateEllipses.end(), [&](const ELLIPSE& candidate){
        if(drawn < numMaxEllipses && existsEllipse(module, candidate)){
            img.draw_ellipse(candidate.a0, candidate.b0, candidate.a, candidate.b, candidate.atanK, color, 1,1);
            ++drawn;
        }
    });
    
    while (!dispSpatial.is_closed() && !acc1Spatial.is_closed())
    {
        
        dispSpatial.wait(dispSpatial,acc1Spatial);
        dispSpatial.display(img);
        // When clicking on the image window.
        if (dispSpatial.button())
        {
            CImg<> win_param(1,3);
            win_param(0,0) = dispSpatial.mouse_x();
            win_param(0,1) = dispSpatial.mouse_y();
            
            dispSpatial.display(img);
            //resultSpatial.display(result);
        }
        
    }
    return 0;

}
