////////////////////////////////////////////////////////////////////////////////
//                          PROJET				                              //
////////////////////////////////////////////////////////////////////////////////

#include "CImg.h"
#include <iostream>
#include <math.h>
#include <vector>
#define DETECT 20

using namespace cimg_library;
using namespace std;
#define DEGREES_TO_RADIANS(degrees) (degrees*M_PI/180)
#define RADIANS_TO_DEGREES(radians) (radians*180/M_PI)

const unsigned char col1[3]={255,255,0}, col2[3]={0,0,0}, col3[3]={255,0,0}, col4[3]={0,255,0};
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

int maxHisto(vector<unsigned long> &v){
    unsigned long max = v.at(0);
    int index = 0;
    for(int i = 0; i < v.size(); ++i){
        if(v.at(i) > max){
            max = v.at(i);
            index = i;
        }
    }
    return index;
    
}
void MaxAccKN(CImg<> &Acc, int& atanK, int& atanN){
    
    int max = -1;
    for(int x = 0; x < Acc.width(); ++x)
        for(int y = 0; y < Acc.height(); ++y){
            if(Acc(x,y) > max){
                max = Acc(x,y);
                atanK = x;
                atanN = y;
            }
        }
}
/*
void checkAngle(float& angle){
    if(angle > 2*M_PI){
        angle = std::fmod(angle, 2*M_PI);
    }
    else if(angle < 0){
        if(angle < -2*M_PI)
            angle = std::fmod(angle, 2*M_PI);
        angle += 2*M_PI;
    }
    
}
*/
/*
void checkAngleAtan(float& x, float& y, float& angle){
    float angleD = RADIANS_TO_DEGREES(angle);
    if(angleD < 0){
        if(x < 0)
            angleD = -angleD;
        else
            angleD += 180;
        angle = DEGREES_TO_RADIANS(angleD);
            
    }
        
}
 */
void checkAngle(float& angle){
   
    if(angle < 0){
        if(angle < -2*M_PI)
            angle = fmod(angle, 2*M_PI) +2*M_PI;
        else
            angle += 2*M_PI;
    }
    else if(angle > 2*M_PI){
        angle = fmod(angle, 2*M_PI);
    }
    
    
    
}

void DrawEllipse(CImg<> &ImgIn, CImg<> &Acc1, CImg<> &Acc2, vector<unsigned long> &Histo){
    
    int a0 = 200;
    int b0 = 200;
    int a,b;
    int ax = maxHisto(Histo);
    float ay, by, bx;
    int atanK, atanN;
    float K, N,atanKrad, atanNrad;
    MaxAccKN(Acc2, atanK, atanN);
    atanKrad = DEGREES_TO_RADIANS(atanK);
    atanNrad = DEGREES_TO_RADIANS(atanN);
    K = tan(atanKrad);
    N = tan(atanNrad);
    ay = K*ax;
    by = N*ax;
    bx = -ay*by /ax;
    a = sqrt(ax*ax + ay*ay);
    b = sqrt(bx*bx+by*by);
    
    
    unsigned char purple[] = { 233,100,125 };
    ImgIn.draw_ellipse(a0, b0, a, b,atanK,purple,1);
}

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


/*
 void EllipseAccumulator(CImg<> ImgIn, CImg<> &Acc1, CImg<> &Acc2, vector<unsigned long> Histo)
 {
 
 Acc1.fill(0);
 Acc2.fill(0);
 
 CImg<> module = Module(ImgIn);
 CImgList<> grad = ImgIn.get_gradient("xy", 3);
 
 cimg_forXY(module, x1, y1)
 {
 
 
 cimg_forXY(module, x2, y2)
 {
 
 
 //if(y2-y1 + x2-x1 < 25)
 //  continue;
 if(sqrt(x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) > 25){
 
 // Initialization
 float xv = -grad.at(1)(x1, y1);
 float yv = grad.at(0)(x1, y1);
 float xw = -grad.at(1)(x2, y2);
 float yw = grad.at(0)(x2, y2);
 
 int Y = y2 - y1;
 int X = x2 - x1;
 // float m1 = -xv/yv;
 //float m2 = -xw/yw;
 float slope1 = -xv/yv;
 float slope2 = -xw/yw;
 float theta1 = tan(slope1);
 float theta2 = tan(slope2);
 float m1 = tan(theta1);
 float m2 = tan(theta2);
 float M1 = m1 * m2;
 float M2 = m1 + m2;
 
 if (abs(m1-m2) > DEGREES_TO_RADIANS(10)) {
 
 //int xp, yp;
 int xm = (x1+x2)/2, ym = (y1+y2)/2;
 
 float xt = ((xv*xw)/(yv*xw - yw*xv)) * (y2-y1 + x1*(yv/xv) - x2*(yw/xw));
 float yt = y1 + (xt - x1) * (yv/xv);
 
 if (xt < 0 || xt >= module.width() || yt < 0 || yt >= module.height())
 continue;
 
 // Calcul de P
 
 for (xp = xm; xp < xt; xp++) {
 yp = ym + (xp-xm)*(yt-ym)/(xt-xm);
 //std::cout << xp << " " << yp << std::endl;
 if(module(xp,yp) > DETECT)
 break;
 }
 
 int b0 = 0;
 
 // Variation de a0
 for (int a0 = 0; a0 < module.width(); ++a0)
 {
 
 // Acc1 value
 //b0 = ( ( a0 - xp ) * ( (2*M1*X) - (Y*M2) ) / (X*M2 - 2*Y) ) + yp;
 b0 = (a0 - xm) * ((yt-ym)/(xt-xm)) + ym;
 
 if (b0 >= module.height() || b0 < 0)
 continue;
 
 Acc1(a0, b0) += 1;
 
 }
 //FIN CENTRES
 
 
 
 
 
 }
 
 }
 
 }
 }
 
 
 }
 */

/******************************************
 Calcul de la phase du gradient de l'image imageIn
 ******************************************/
CImg<float> Phase (CImg<unsigned char> imageIn)
{
    CImg <> Px(3,3,1,1), Py(3,3,1,1);
    
    Px(0,0)=-1; Px(0,1)=-1; Px(0,2)=-1;
    Px(1,0)=0;  Px(1,1)=0;  Px(1,2)=0;
    Px(2,1)=1;  Px(2,0)=1;  Px(2,2)=1;
    
    Py(0,0)=-1; Py(0,1)=0; Py(0,2)=1;
    Py(1,0)=-1; Py(1,1)=0; Py(1,2)=1;
    Py(2,1)=-1; Py(2,0)=0; Py(2,2)=1;
    
    CImg<double> img_x = imageIn.get_convolve(Px);
    CImg<double> res   = imageIn.get_convolve(Py);
    res.atan2(img_x);
 
    
    return res;
}

void EllipseAccumulator(CImg<> ImgIn, CImg<> &Acc1, CImg<> &Acc2, vector<unsigned long> Histo)
{
    
	Acc1.fill(0);
	Acc2.fill(0);
    
	CImg<> module = Module(ImgIn);
	CImgList<> grad = ImgIn.get_gradient("xy", 3);
    
	cimg_forXY(module, x1, y1)
	{
        
        cimg_forXY(module, x2, y2)
        {
            if(y2-y1 + x2-x1 < 25)
                continue;
            
            // Initialization
            float xv = -grad.at(1)(x1, y1);
            float yv = grad.at(0)(x1, y1);
            float xw = -grad.at(1)(x2, y2);
            float yw = grad.at(0)(x2, y2);
            
            float m1 = -xv/yv;
            float m2 = -xw/yw;
            
            if ( module(x1, y1) > DETECT && module(x2, y2) > DETECT && abs(m1-m2) > 10) {
                
                int xp, yp;
                int xm = (x1+x2)/2, ym = (y1+y2)/2;
                
                float xt = ((xv*xw)/(yv*xw - yw*xv)) * (y2-y1 + x1*(yv/xv) - x2*(yw/xw));
                float yt = y1 + (xt - x1) * (yv/xv);
                
                if (xt < 0 || xt >= module.width() || yt < 0 || yt >= module.height())
                    continue;
                
                
                // Calcul de P
                for (xp = xm; xp < xt; xp++) {
                    yp = ym + (xp-xm)*(yt-ym)/(xt-xm);
                    //std::cout << xp << " " << yp << std::endl;
                    if(module(xp,yp) > DETECT)
                        break;
                }
                
                int b0 = 0;
                
                // Variation de a0
                for (int a0 = 0; a0 < module.width(); ++a0)
                {
                    
                    // Acc1 value
                    //b0 = ( ( a0 - xp ) * ( (2*M1*X) - (Y*M2) ) / (X*M2 - 2*Y) ) + yp;
                    b0 = (a0 - xm) * ((yt-ym)/(xt-xm)) + ym;
                    
                    if (b0 >= module.height() || b0 < 0)
                        continue;
                    
                    Acc1(a0, b0) += 1;
                    
                }
                
            }
            
        }
        
    }
    
}




void AccumulateKN(CImg<> &ImgIn, CImg<> &Acc1, CImg<> &Acc2, vector<unsigned long>& Histo, int minA, int maxA, int minB, int maxB){
    int
    xp, yp, atanK,atanN, a, b;
    
    float
    slopeFirstDerivativeP, slopeSecondDerivativeP, angleDerivativeP1, angleDerivativeP2, angleFirstDerivativeP, angleSecondDerivativeP, m1, m2, M1, M2, N, K, x0, y0,phi1, phi2, ro, X, Y,ax, ay, bx, by;
    float margin = 0.0001;
    CImg<> tangent(ImgIn.width(), ImgIn.height());
    CImg<> module(ImgIn.width(), ImgIn.height());
    module = Module(ImgIn);
    CImg<> phase(ImgIn.width(), ImgIn.height());
    phase = Phase(ImgIn);
    CImgList<> grad = ImgIn.get_gradient("xy", 3);
    
   // CImgList<> grad = ImgIn.get_gradient("xy", 3);
    int a0 = 200;
    int b0 = 200;

    cimg_forXY(module, x1, y1)
        if(module(x1, y1) > 0) //There is an edge point
            cimg_forXY(module, x2, y2)
            if(module(x2, y2) > 0 && sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)) > 25){
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
                    xp = ((-mt)*(xt) - b0 + yt + slopeSecondDerivativeP*a0)/(slopeSecondDerivativeP - mt);
                    yp = slopeSecondDerivativeP*(xp - a0) + b0;
                    //for(xp = 0; xp < ImgIn.width(); ++xp){
                      //  yp = floor(slopeSecondDerivativeP*(xp - a0) + b0 +0.5f);
                        if(xp >0 && xp < ImgIn.width() && yp > 0 && yp < ImgIn.height() && module(xp, yp) > 0){
                            phi1 = atan2(Y,X);
                            //phi2 = atan2(2*M1*X - Y*M2, X*M2-2*Y);
                            phi2 = atan((2*M1*X - Y*M2)/(X*M2-2*Y));
                          
                          
                           // ro = phase(xp, yp);
                            for(int i = 0; i < 180; i+=5){
                                ro = DEGREES_TO_RADIANS(i);
                            
                         //   if(abs(ro-90) <= margin || abs(ro-DEGREES_TO_RADIANS(270)) <= margin) continue;
                                if(i == 0 || i == 90 ) continue;
                            if(abs((phi1-ro)-DEGREES_TO_RADIANS(90)) <= margin || abs((phi1-ro)-DEGREES_TO_RADIANS(270)) <= margin || abs((phi2-ro) - DEGREES_TO_RADIANS(90)) <= margin || abs((phi2-ro)-DEGREES_TO_RADIANS(270)) <= margin) continue;
                                N = sqrt(abs(tan(phi1 - ro)*tan(phi2 - ro)));
                                K = tan(ro);
                                x0 = (xp - a0)/sqrt(K*K+1) + ((yp - b0)*K)/sqrt(K*K+1);
                                y0 = ((xp - a0)*K)/sqrt(K*K+1) + (yp - b0)/sqrt(K*K+1);
                                if((x0 + a0) > 0 && (x0 + a0) < ImgIn.width() && (y0 + b0) > 0 && (y0 +b0) < ImgIn.height() && module(x0+a0, y0 + b0) > 0){
                                ax = sqrt(abs((y0*y0 + x0*x0*N*N)/(N*N*(1 + K*K))));
                                if(ax <= margin) continue;
                                        ay = K*ax;
                                        by = N*ax;
                                        bx = -ay*by/ax;
                                        a = sqrt(ax*ax + ay*ay);
                                        b = sqrt(bx*bx + by*by);
                                        if(a > minA && a < maxA && b > minB && b < maxB){
                                            atanN = RADIANS_TO_DEGREES(atan(N));
                                            if(atanN < 0)
                                                atanN +=360;
                                            
                                            
                                            Acc2(i, atanN) += 1;
                                            std::cout << "Added to Acc2:" << i << ", " << atanN << std::endl;
                                            ++Histo[floor(ax+0.5f)];
                                            std::cout << "Added ax: " << ax << std::endl;
                                            
                                        }
                                }
                        }
                        
                        }
                            }
                        }
}






    


    
  


/*******************************************************************************
 
 Main
 
 *******************************************************************************/
int main(int argc,char **argv)
{
    cimg_usage("Retrieve command line arguments");
    const char* filename = cimg_option("-i","/Users/rubcuevas/Desktop/Algorithmie de l'image/EllipseDetection/EllipseMerged/ellipsefilled.bmp","Input image file");
    
    const int minA = 12;
    const int maxA = 17;
    const int minB = 28;
    const int maxB = 32;
    // Opening of filename
    CImg<> img(filename);
    
    //CImg<> img(400, 400, 1, 3, 255);
    //img.draw_ellipse(200, 200, 50, 100, cimg::PI/4, col1);
    
    // Accumulator Computation
    CImg<> Acc1(img.width(), img.height());
    CImg<> Acc2(180, 180);
    std::vector<unsigned long> Histo(img.width()/2);
    
    //EllipseAccumulator(img, Acc1, Acc2, Histo);
    unsigned long Acc1Values[Acc1.width()][Acc1.height()];
    for(int i = 0; i < Acc1.width(); ++i)
        for(int j = 0; j < Acc1.height(); ++j)
            Acc1Values[i][j] = Acc1(i,j);
    AccumulateKN(img, Acc1, Acc2, Histo, minA, maxA, minB, maxB);
    
    // Display
    CImg<> result(img.width(), img.height());
    CImgDisplay dispSpatial(img,"Input Image");
    CImgDisplay acc1Spatial(Acc1,"Accumulator 1 map");
    CImgDisplay acc2Spatial(Acc2,"Accumulator 2 map");
    CImgDisplay resultSpatial(result, "Result");
    int maxH = maxHisto(Histo);
    CImg<> Acc2Thresholded(Acc2.width(), Acc2.height());
    //MaxDetection(Acc2, Acc2Thresholded);
    DrawEllipse(result, Acc1, Acc2, Histo);
    
    while (!dispSpatial.is_closed() && !acc1Spatial.is_closed() && !acc2Spatial.is_closed())
    {
        
        dispSpatial.wait(dispSpatial,acc1Spatial);
        
        // When clicking on the image window.
        if (dispSpatial.button())
        {
            CImg<> win_param(1,3);
            win_param(0,0) = dispSpatial.mouse_x();
            win_param(0,1) = dispSpatial.mouse_y();
            
            dispSpatial.display(img);
            resultSpatial.display(result);
        }
        
    }
    return 0;
}

