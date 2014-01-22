////////////////////////////////////////////////////////////////////////////////
//                          PROJET				                              //
////////////////////////////////////////////////////////////////////////////////

#include "CImg.h"
#include <iostream>
#include <math.h>
#define DETECT 50

using namespace cimg_library;

const unsigned char col1[3]={255,255,0}, col2[3]={0,0,0}, col3[3]={255,0,0}, col4[3]={0,255,0};

CImg<> Module (CImg<> image)
{
	CImgList<float> list = image.get_gradient("xy", 0);
	return (list.at(0).sqr()+list.at(1).sqr()).get_sqrt();
}

void EllipseAccumulator(CImg<> ImgIn, CImg<> &Acc1, CImg<> &Acc2, int* Histo)
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
                
                /*for (int K = 0; K < 360; ++K) {
                    
                    // Acc2 value
                    float phi1 = atan2(Y, X);
                    float phi2 = atan2( (2*M1*X) - (Y*M2) , (X*M2 - 2*Y) );
                    //float rho = atan2(K*cimg::PI, 360);
                    float rho = K;
                    float N2 = tan(phi1 - rho) * tan(phi2 - rho);
                    float N = atan(sqrt(N2));
                    
                    if (N < 0)
                        N += cimg::PI;
                    
                    //std::cout << N << std::endl;
                    
                    Acc2(K, N) += 1;
                    
                    // Histo value
                    //float x0 = ( (xp - a0) / (sqrt(K*K+1)) ) + ( (yp - b0)*K / (sqrt(K*K+1)) );
                    //float y0 = ( (xp - a0)*K / (sqrt(K*K+1)) ) + ( (yp - b0) / (sqrt(K*K+1)) );
                    //int ax = sqrt( (y0*y0 + (x0*x0+N2)) / (N2*(1+K*K)) );
                    
                    //Histo[ax] += 1;
                    
                    
                }*/
                
                
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
    const char* filename = cimg_option("-i","/Users/rubcuevas/Desktop/Algorithmie de l'image/EllipseDetection/EllipseDetection/ellipse.bmp","Input image file");
    
    // Opening of filename
    CImg<> img(filename);
    
    //CImg<> img(400, 400, 1, 3, 255);
    //img.draw_ellipse(200, 200, 50, 100, cimg::PI/4, col1);
    
    // Accumulator Computation
    CImg<> Acc1(img.width(), img.height());
    CImg<> Acc2(360, 360);
    int* Histo = new int[img.width()];
    
    EllipseAccumulator(img, Acc1, Acc2, Histo);
    
    // Display
    CImgDisplay dispSpatial(img,"Input Image");
    CImgDisplay acc1Spatial(Acc1,"Accumulator 1 map");
    CImgDisplay acc2Spatial(Acc2,"Accumulator 2 map");
    
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
