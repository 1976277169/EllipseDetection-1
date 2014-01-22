////////////////////////////////////////////////////////////////////////////////
//                          PROJET				                              //
////////////////////////////////////////////////////////////////////////////////

#include "CImg.h"
#include <iostream>
#include <math.h>
#include <vector>
#define DETECT 50

using namespace cimg_library;
using namespace std;
#define DEGREES_TO_RADIANS(degrees) (degrees*M_PI/180)
#define RADIANS_TO_DEGREES(radians) (radians*180/M_PI)

const unsigned char col1[3]={255,255,0}, col2[3]={0,0,0}, col3[3]={255,0,0}, col4[3]={0,255,0};

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
    cimg_forXY(Acc, x, y){
        if(Acc(x,y) > max){
            max = Acc(x,y);
            atanK = x;
            atanN = y;
        }
    }
}

void DrawEllipse(CImg<> &ImgIn, CImg<> &Acc1, CImg<> &Acc2, vector<unsigned long> &Histo){
    int a0 = 200;
    int b0 = 200;
    int a,b;
    int ax = maxHisto(Histo);
    int ay, by, bx;
    int atanK, atanN, K, N;
    MaxAccKN(Acc2, atanK, atanN);
    K = tan(atanK);
    N = tan(atanN);
    ay = K*ax;
    by = N*ax;
    bx = -ay*-by /ax;
    a = sqrt(ax*ax + ay*ay);
    b = sqrt(bx*bx+by*by);
   
    unsigned char purple[] = { 233,0,125 };
    ImgIn.draw_ellipse(a0, b0, a, b, DEGREES_TO_RADIANS(atanK),purple,1);
}

CImg<> Module (CImg<> image)
{
	//CImgList<float> list = image.get_gradient("xy", 0);
    CImgList<float> G = image.get_gradient("xy", 3);
    CImg<> module = (G[0].sqr()+G[1].sqr()).get_sqrt();
    CImg<> ToReturn(image.width(), image.height());
    cimg_forXY(image, x, y){
        if(module(x,y) > DETECT)
            ToReturn(x,y) = module(x,y);
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
            
            int Y = y2 - y1;
            int X = x2 - x1;
            float m1 = -xv/yv;
            float m2 = -xw/yw;
            float M1 = m1 * m2;
            float M2 = m1 + m2;
            
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
void AccumulateKN(CImg<> &ImgIn, CImg<> &Acc1, CImg<> &Acc2, vector<unsigned long>& Histo){
     CImg<> module(ImgIn.width(), ImgIn.height());
    module = Module(ImgIn);
    CImgList<> grad = ImgIn.get_gradient("xy", 3);
    //CImg<> tangent = ImgIn.get_tan();
    //CImg<> tangent = module.get_tan();
    //cimg_forXY(module, x1, y1){
    int a0 = 200;
    int b0 = 200;
    for(int x1 = 0; x1 < module.width(); ++x1){
        for(int y1 = 0; y1 < module.height(); ++y1){
          if(module(x1,y1) > 0){
              
              
              for(int x2 = 0; x2 < module.width(); ++x2){
                  for(int y2 = 0; y2 < module.height(); ++y2)
            {
                
                if(module(x2,y2) > 0){
                    ////
                    if(sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)) > 25){
                        
                        float X = x2 - x1;
                        float Y = y2 - y1;
                        float angleFirstDerivativeP1 = atan2(grad[1](x1,y1),grad[0](x1,y1));
                        if(angleFirstDerivativeP1 < 0)
                            angleFirstDerivativeP1 = -angleFirstDerivativeP1;
                        float angleFirstDerivativeP2 = atan2(grad[1](x2,y2),grad[0](x2,y2));
                        if(angleFirstDerivativeP2 < 0)
                            angleFirstDerivativeP2 = -angleFirstDerivativeP2;
                        float m1 = tan(angleFirstDerivativeP1);
                        float m2 = tan(angleFirstDerivativeP2);
                        if(abs(m1 -m2) < DEGREES_TO_RADIANS(10))
                            continue;
                        
                        float M1 = m1 * m2;
                        float M2 = m1 + m2;
                        float angleFirstDerivativeP = Y/X;
                        checkAngle(angleFirstDerivativeP);
                        float angleSecondDerivativeP = ((2*M1*X) - (Y*M2)) / ((X*M2 - 2*Y));
                        checkAngle(angleSecondDerivativeP);
                        
                        float phi1 = atan(angleFirstDerivativeP);
                        checkAngle(phi1);
                        float phi2 = atan(angleSecondDerivativeP);
                        checkAngle(phi2);
                        float ro;
                        float N, K, x0, y0, ax;
                        int xp, yp;
                        for(xp = a0;xp < module.width(); ++xp){
                             yp = angleSecondDerivativeP * (xp - a0) + b0;
                            if(yp < module.height()){
                                
                               
                            if(module(xp, yp) > 0){
                            for(int i = 0; i < 180; ++i){
                                ro = DEGREES_TO_RADIANS(i);
                                checkAngle(ro);
                                K = tan(ro);
                                if(K < 0) K = -K;
                                
                                N = sqrt(abs(tan(phi1 - ro)*tan(phi2 - ro)));
                                x0 = ((xp - a0)/(sqrt(K*K + 1)) + ((yp - b0)*K)/sqrt(K*K+1));
                                y0 = ((xp - a0)*K)/(sqrt(K*K+1))+((yp - b0)/(sqrt(K*K + 1)));
                                ax = sqrt(abs((y0*y0 + x0*x0*N*N)/(N*N*(1+K*K))));
                                int atanN = atan(N);
                                checkAngle(N);
                                Acc2(i, atanN);
                                ++Histo[ax];
                        }
                            }
                        }
                        }
                    }
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
    const char* filename = cimg_option("-i","/Users/rubcuevas/Desktop/Algorithmie de l'image/EllipseDetection/EllipseMerged/ellipse.bmp","Input image file");
    
    // Opening of filename
    CImg<> img(filename);
    
    //CImg<> img(400, 400, 1, 3, 255);
    //img.draw_ellipse(200, 200, 50, 100, cimg::PI/4, col1);
    
    // Accumulator Computation
    CImg<> Acc1(img.width(), img.height());
    CImg<> Acc2(360, 360);
    std::vector<unsigned long> Histo(img.width()/2);
    
    //EllipseAccumulator(img, Acc1, Acc2, Histo);
    unsigned long Acc1Values[Acc1.width()][Acc1.height()];
    for(int i = 0; i < Acc1.width(); ++i)
        for(int j = 0; j < Acc1.height(); ++j)
            Acc1Values[i][j] = Acc1(i,j);
    AccumulateKN(img, Acc1, Acc2, Histo);
    
    // Display
    CImg<> result(img.width(), img.height());
    CImgDisplay dispSpatial(img,"Input Image");
    CImgDisplay acc1Spatial(Acc1,"Accumulator 1 map");
    CImgDisplay acc2Spatial(Acc2,"Accumulator 2 map");
    CImgDisplay resultSpatial(result, "Result");
    int maxH = maxHisto(Histo);
    DrawEllipse(img, Acc1, Acc2, Histo);
    
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
        }
        
    }
    return 0;
}

