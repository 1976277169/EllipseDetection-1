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
void EllipseAccumulator(CImg<> ImgIn, CImg<> &Acc1, CImg<> &Acc2, vector<unsigned long> Histo)
{

	Acc1.fill(0);
	Acc2.fill(0);

	CImg<> module = Module(ImgIn);
	CImgList<> grad = ImgIn.get_gradient("xy", 3);

	cimg_forXY(module, x1, y1)
	{
        std::cout <<y1 << std::endl;
	
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
                
                int xp, yp;
                int xm = (x1+x2)/2, ym = (y1+y2)/2;
                
                float xt = ((xv*xw)/(yv*xw - yw*xv)) * (y2-y1 + x1*(yv/xv) - x2*(yw/xw));
                float yt = y1 + (xt - x1) * (yv/xv);
                
                if (xt < 0 || xt >= module.width() || yt < 0 || yt >= module.height())
                    continue;
                
                // Calcul de P
                /*
                for (xp = xm; xp < xt; xp++) {
                    yp = ym + (xp-xm)*(yt-ym)/(xt-xm);
                    //std::cout << xp << " " << yp << std::endl;
                    if(module(xp,yp) > DETECT)
                        break;
                }
                */
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

void AccumulateKN(CImg<> &ImgIn, CImg<> &Acc1, CImg<> &Acc2, vector<unsigned long>& Histo){
     CImg<> module(ImgIn.width(), ImgIn.height());
    module = Module(ImgIn);
    CImgList<> grad = ImgIn.get_gradient("xy", 3);
    //CImg<> tangent = ImgIn.get_tan();
    //CImg<> tangent = module.get_tan();
    //cimg_forXY(module, x1, y1){
    for(int x1 = 0; x1 < module.width(); ++x1){
        for(int y1 = 0; y1 < module.height(); ++y1){
          if(module(x1,y1) > 0){
            //float theta1 = tangent(x1,y1);
            //float slope1 = tan(theta1);
           // cimg_forXY(module,x2,y2)
              
              
              for(int x2 = 0; x2 < module.width(); ++x2){
                  for(int y2 = 0; y2 < module.height(); ++y2)
            {
                
                if(module(x2,y2) > 0){
                    ////
                    if(sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)) > 25){
                    float xv = -grad[1](x1, y1);
                    float yv = grad[0](x1, y1);
                    float xw = -grad[1](x2, y2);
                    float yw = grad[0](x2, y2);
                    float slope1 = -xv/yv;
                    float slope2 = -xw/yw;
                    float theta1 = tan(slope1);
                    float theta2 = tan(slope2);
                    float t1 = ((xv*xw)/(yv*xw - yw*xv)) * (y2-y1 + x1*(yv/xv) - x2*(yw/xw));
                    float t2 = y1 + (t1 - x1) * (yv/xv);
                    
                ////
                //float theta2 = tangent(x2,y2);
                //float slope2 = tan(theta2);
        if(abs((tan(theta1) - tan(theta2))) > DEGREES_TO_RADIANS(10)){
    
    //float t1 = (y1 - y2 - x1*slope1 + x2*slope2) / (slope2-slope1);
    //float t2 = (slope1*slope2*(x2-x1) - y2*slope1 + y1*slope2) / (slope2-slope1);
    float m1 = (x1 + x2) / 2;
    float m2 = (y1 + y2) / 2;
    float slope = (t2-m2) / (t1-m1);
    float b = (m2*t1 - m1*t2) / (t1-m1);
    float xp, yp; // Coordinates of P
    
    float dp, d2p, M1, M2, X, Y, K, N, phi1, phi2, x0, y0;
    int ax;
    //cimg_forXY(Acc1, a0, b0){
        //Revise condition!
     //   if(Acc1(a0,b0) > 100){
            //Calculate point P
            int a0 = 200;
            int b0 = 200;
            X = x2 - x1;
            Y = y2 - y1;
            dp = Y/X;
            //M1 = slope1 * slope2;
            //M2 = slope1 + slope2;
            M1 = tan(theta1) * tan(theta2);
            M2 = tan(theta1) + tan(theta2);
            d2p = (2*M1*X - Y*M2)/(X*M2 - 2*Y);
           //Calculer p
            
            float anglep1 = atan2(y1,x1);
            float anglep2 = atan2(y2,x2);
            float anglepp = (anglep1 + anglep2)/2;
            
            
            
           // xp = (d2p*a0 - b0) / (d2p - tan((theta1 + theta2) / 2));
           // yp = d2p*xp -d2p*a0 + b0;
            xp = (d2p*a0 - b0)/(d2p - tan(anglepp));
            yp = xp*tan(anglepp);
                if(xp >= 0 && xp < ImgIn.width() && yp >= 0 && yp <= ImgIn.height() && module(xp, yp) > 0   ){
                    //Calculate parameters
                    phi1 = atan2(Y,X);
                    phi2 = atan2((2*M1*X - Y*M2),(X*M2 - 2*Y));
                    if(phi1 < 0)
                        phi1 = -phi1;
                    if(phi2 < 0)
                        phi2 = -phi2;
                    //Generate ro
                    for(int ro = 0; ro < 180; ro++){
                        //if(ro == 90){
                          //  ro += 10;
                            //continue;
                    //}
                        float radians = DEGREES_TO_RADIANS(ro);
                        N = sqrt(abs(tan(phi1 - radians) - tan(phi2 - radians)));
                        K = tan(radians);
                        x0 = ((xp - a0)/(sqrt(K*K + 1))) + ((yp - b0)*K)/sqrt(K*K + 1);
                        y0 = (((xp - a0) * K) /(((sqrt(K*K + 1))) + (yp - b0))/sqrt(K*K + 1));
                        if(xp >= 0 && xp < ImgIn.width() && yp >= 0 && yp <= ImgIn.height()){
                            float divi = sqrt((y0*y0 + x0*x0*N*N));
                            float divis = sqrt(((N*N)*(1+K*K)));
                            ax = floor(divi/divis + 0.5f);
                            int atanN = atan(N);
                            if(atanN < 0)
                                atanN = atanN + 2*M_PI;
                            atanN = RADIANS_TO_DEGREES(atanN);
                            
                           // int atanN = floor(RADIANS_TO_DEGREES(atan(N)) + 0.5f);
                            //int atanK = floor(RADIANS_TO_DEGREES(atan(K)) + 0.5f);
                           
                            if(ax > 0 && ax < ImgIn.width()/2){
                                    Acc2(ro, atanN)  += 1;
                            
                            
                                    ++Histo[ax];
                                
                                ImgIn.draw_point(xp, yp, col1, 0.01);
                            }
                        }
                    }
                
            }
            }
                }
            //}
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
    CImg<> Acc2(180, 360);
    std::vector<unsigned long> Histo(img.width());
    
    //EllipseAccumulator(img, Acc1, Acc2, Histo);
    AccumulateKN(img, Acc1, Acc2, Histo);
    
    // Display
    CImgDisplay dispSpatial(img,"Input Image");
    CImgDisplay acc1Spatial(Acc1,"Accumulator 1 map");
    CImgDisplay acc2Spatial(Acc2,"Accumulator 2 map");
    int maxH = maxHisto(Histo);
    
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

