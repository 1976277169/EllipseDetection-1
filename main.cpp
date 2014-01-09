////////////////////////////////////////////////////////////////////////////////
//                          PROJET				                              //
////////////////////////////////////////////////////////////////////////////////

#include "CImg.h"
#include <math.h>
#include <vector>

#define DETECT 50
#define DEGREES_TO_RADIANS(degrees) (degrees*M_PI / 180)
#define RADIANS_TO_DEGREES(radians) (radians/M_PI*180)
using namespace cimg_library;
using namespace std;

/******************************************************************************* 
 ADDED BY ME
 
 DrawHoughLine: Draw the line of the Hough Transform and circle in the accumulator
 
 ImgIn : Input Image where to draw the lines
 AccIn : Accumulator image where to draw circles
 param : Line in parametric form
 
 *******************************************************************************/
void DrawHoughLine(CImg<> &ImgIn,CImg<> &AccIn,CImg<> param,int NbTheta, int NbRho)
{
    const unsigned char col1[3]={255,255,0}, col2[3]={0,0,0};
    const float col3[1]={AccIn.max() ? AccIn.max() : 1 };
    
    cimg_forX(param,i)
    {
        const double
        rhomax   = sqrt(ImgIn.width()*ImgIn.width()+ImgIn.height()*ImgIn.height())/2.,
        thetamax = 2.*cimg::PI,
        rho   = param(i,1)*rhomax/NbRho,
        theta = param(i,0)*thetamax/NbTheta,
        x = ImgIn.width()/2 + rho*std::cos(theta),
        y = ImgIn.height()/2 + rho*std::sin(theta);
        
        const int
        x0 = (int)(x+1000*std::sin(theta)),
        y0 = (int)(y-1000*std::cos(theta)),
        x1 = (int)(x-1000*std::sin(theta)),
        y1 = (int)(y+1000*std::cos(theta));
        
        ImgIn.
        draw_line(x0,y0,x1,y1,col1,1.0f,0xF0F0F0F0).draw_line(x0,y0,x1,y1,col2,1.0f,0x0F0F0F0F).
        draw_line(x0+1,y0,x1+1,y1,col1,1.0f,0xF0F0F0F0).draw_line(x0+1,y0,x1+1,y1,col2,1.0f,0x0F0F0F0F).
        draw_line(x0,y0+1,x1,y1+1,col1,1.0f,0xF0F0F0F0).draw_line(x0,y0+1,x1,y1+1,col2,1.0f,0x0F0F0F0F);
        
        AccIn.
        draw_circle(param(i,0),param(i,1),7,col3,1.0f,0xF0F0F0F0).
        draw_circle(param(i,0),param(i,1),8,col3,1.0f,0xF0F0F0F0);
    }
}

void Module(CImg<> ImgIn, CImg<> &ImgOut){
    CImgList<> grad = ImgIn.get_gradient("xy", 3);
    cimg_forXY(ImgIn, x, y){
        double mod = sqrt(grad[0](x,y)*grad[0](x,y) + grad[1](x,y)*grad[1](x,y));
            if(mod > DETECT){
                ImgOut(x,y) = mod;
            }
            else{
                ImgOut(x,y) = 0;
            }
    }
    
    
}

void AccumulateCenter(CImg<> ImgIn, CImg<> &Acc1, float x1, float x2, float y1, float y2, float slope1, float slope2){
    /*
     * Intersection point of tangents: (t1, t2).
     * Midpoint between (x1,y1) and (x2,y2) : (m1, m2)
     */
    float t1 = (y1 - y2 - x1*slope1 + x2*slope2) / (slope2-slope1);
    float t2 = (slope1*slope2*(x2-x1) - y2*slope1 + y1*slope2) / (slope2-slope1);
    float m1 = (x1 + x2) / 2;
    float m2 = (y1 + y2) / 2;
    
    float slope = (t2-m2) / (t1-m1);
    float b = (m2*t1 - m1*t2) / (t1-m1);
    int x_h, x_end, y_h, y_end;
    
    //////////
    if (abs(t2-m2) < abs(t1-m1)){
    /* If (m1,m2) is to the left of (t1,t2)
     */
    if (m1 < t1)
    {
        if (slope == 0)
        {
            x_end = 0;
        }
        else if (slope > 0)
        {
            /* zero intercept
             * i.e. y = 0
             */
            x_end = (int )(-b/slope);
        }
        else /* if slope < 0 */
        {
            /* y = edge of image
             * i.e. y = aImage->y
             */
            x_end = (int )( (((double )ImgIn.height())-b) / slope );
        }
        
        
        /* Increase x_end by one.
         * Otherwise, because of rounding errors,
         * we may get a value of y_h=aImage->y
         * which would give a go outside the bounds
         * of the array histogram[][].
         */
        x_end++;
        
        
        if (x_end < 0)
            x_end = 0;
        
        
        for (x_h=(int )m1; x_h > x_end; x_h--)
        {
            y_h = (int )(slope*(double )x_h + b);
            Acc1(x_h,y_h)++;
        }
    }
    else    /* if m1 > t1 */
    {
        /* If (m1,m2) is to the right of (t1,t2)
         */
        if (slope == 0)
        {
            x_end = ImgIn.width();
        }
        else if (slope > 0)
        {
            /* y = edge of image
             * i.e. y = aImage->y
             */
            x_end = (int )( (((double )ImgIn.height())-b) / slope );
        }
        else /* if slope < 0 */
        {
            /* zero intercept
             * i.e. y = 0
             */
            x_end = (int )(-b/slope);
        }
        
        
        /* Decrease x_end by one.
         * Otherwise, because of rounding errors,
         * we may get a value of y_h=aImage->y
         * which would give a go outside the bounds
         * of the array histogram[][].
         */
        x_end--;
        
        
        if (x_end > ImgIn.width())
            x_end = ImgIn.width();
        
        for (x_h=(int )m1; x_h < x_end; x_h++)
        {
            y_h = (int )(slope*(double )x_h + b);
            Acc1(x_h,y_h)++;
        }
        
    }
    }
    else{
        /*
         * y co-ord varies more quickly
         * than x co-ord (i.e. slope > 1),
         * so index the line by the y co-ord.
         */
        slope = (t1-m1) / (t2-m2);
        b = (m1*t2 - m2*t1) / (t2-m2);
        
        /* If (m1,m2) is below (t1,t2)
         */
        if (m2 < t2)
        {
            if (slope == 0)
            {
                y_end = 0;
            }
            else if (slope > 0)
            {
                /* zero intercept
                 * i.e. x = 0
                 */
                y_end = (int )(-b/slope);
            }
            else /* if slope < 0 */
            {
                /* x = edge of image
                 * i.e. x = aImage->x
                 */
                y_end = (int )( (((double )ImgIn.width())-b) / slope );
            }
            
            
            /* Increase y_end by one.
             * Otherwise, because of rounding errors,
             * we may get a value of x_h=aImage->x
             * which would give a go outside the bounds
             * of the array histogram[][].
             */
            y_end++;
            
            
            if (y_end < 0)
                y_end = 0;
            
            for (y_h=(int )m2; y_h > y_end; y_h--)
            {
                x_h = (int )(slope*(double )y_h + b);
                Acc1(x_h,y_h)++;
            }
        }
        else    /* if m2 > t2 */
        {
            /* If (m1,m2) is above (t1,t2)
             */
            if (slope == 0)
            {
                y_end = ImgIn.height();
            }
            else if (slope > 0)
            {
                /* x = edge of image
                 * i.e. x = aImage->x
                 */
                y_end = (int )( (((double )ImgIn.width())-b) / slope );
            }
            else /* if slope < 0 */
            {
                /* zero intercept
                 * i.e. x = 0
                 */
                y_end = (int )(-b/slope);
            }
            
            
            /* Decrease y_end by one.
             * Otherwise, because of rounding errors,
             * we may get a value of x_h=aImage->x
             * which would give a go outside the bounds
             * of the array histogram[][].
             */
            y_end--;
            
            
            if (y_end > ImgIn.height())
                y_end = ImgIn.height();
            
            for (y_h=(int )m2; y_h < y_end; y_h++)
            {
                x_h = (int )(slope*(double )y_h + b);
                Acc1(x_h,y_h)++;
            }
        }
    }
    
}  /* end of if ((theta1-theta2 > M_PI/10.0)... */
   /* end of for loop based on hp2 */




int FindMaxHisto(vector<int> &Histo){
    int max = -1;
    int pos = -1;
    for(int i = 0; i < sizeof(Histo)/sizeof(int); ++i){
        if(Histo[i] > max){
            max = Histo[i];
            pos = i;
        }
    }
    return pos;
    
}

void AccumulateKN(CImg<> ImgIn, CImg<> &Acc1, CImg<> &Acc2, vector<int> &Histo, float x1, float x2, float y1, float y2, float slope1, float slope2){
    float t1 = (y1 - y2 - x1*slope1 + x2*slope2) / (slope2-slope1);
    float t2 = (slope1*slope2*(x2-x1) - y2*slope1 + y1*slope2) / (slope2-slope1);
    float m1 = (x1 + x2) / 2;
    float m2 = (y1 + y2) / 2;
    float slope = (t2-m2) / (t1-m1);
    float b = (m2*t1 - m1*t2) / (t1-m1);
    float xp, yp; // Coordinates of P
    
    float dp, d2p, M1, M2, X, Y, K, N, phi1, phi2, x0, y0;
    int ax;
    cimg_forXY(Acc1, a0, b0){
        //Revise condition!
        if(Acc1(a0,b0) > 0){
            //Calculate point P
            X = x2 - x1;
            Y = y2 - y1;
            dp = Y/X;
            M1 = slope1 * slope2;
            M2 = slope1 + slope2;
            d2p = (2*M1*X - Y*M2)/(X*M2 - 2*Y);
            /*
             * By formula (4): yp = d2p * (x - a0) +b0
             * P is the intersection of the line MT with the ellipse:
             * mx + b = d2p * x -a0 + b0
             *
             */
            if(abs(slope - d2p > 0.8)){
                xp = (-d2p * a0 + b0 - b)/(slope - d2p);
                yp = d2p*(xp-a0) + b0;
                if(xp >= 0 && xp < ImgIn.width() && yp >= 0 && yp <= ImgIn.height()){
                    //Calculate parameters
                    phi1 = atan(dp);
                    phi2 = atan(d2p);
                    //Generate ro
                    for(int ro = 0; ro < 360; ++ro){
                        if(ro == 90){
                            ro += 10;
                            continue;
                        }
                        float radians = DEGREES_TO_RADIANS(ro);
                        N = sqrt(abs(tan(phi1 - radians) - tan(phi2 - radians)));
                        K = tan(radians);
                        x0 = ((xp - a0)/(sqrt(K*K + 1)) + (yp -b0)*K/sqrt(K*K + 1));
                        y0 = ((xp - a0) * K /(sqrt(K*K + 1) + (yp - b0)/(sqrt(K*K + 1))));
                        if(xp >= 0 && xp < ImgIn.width() && yp >= 0 && yp <= ImgIn.height()){
                            ax = (int)sqrt((y0*y0 + x0*x0*N*N) / (N*N*(1+K*K)));
                            int atanN = RADIANS_TO_DEGREES(floor(atan(N) + 0.5f));
                            int atanK = RADIANS_TO_DEGREES(floor(atan(K) + 0.5f));
                            if(atanK < 0)
                                atanK = atanK + 360;
                            Acc2(atanN, atanK)  += 1;
                            ++Histo[ax];
                        }
                    }
                
            }
            }
                
        }
    }
    
}

void EllipseAccumulator(CImg<> ImgIn, CImg<> &Acc1, CImg<> &Acc2, vector<int> &Histo, float alpha, float sigma)
{
	Acc1.fill(0);
	Acc2.fill(0);
	
    CImg<> module(ImgIn.width(), ImgIn.height());
    Module(ImgIn, module);
    CImg<> tangent = ImgIn.get_tan();
    cimg_forXY(module, x1, y1){
        //There is an edge
        if(module(x1,y1) > 0){
            float theta1 = tangent(x1,y1);
            float slope1 = tan(theta1);
            cimg_forXY(module,x2,y2)
            {
                float theta2 = tangent(x2,y2);
                float slope2 = tan(theta2);
                //Revise////
                if(theta2 > theta1)
                    theta2 -= M_PI;
                /////////
                //Revise condition
                if(abs((theta1 - theta2)) > 10){
                    //Center
                    AccumulateCenter(module, Acc1, x1,x2,y1,y2,slope1,slope2);
                    //Accumulate KN
                    AccumulateKN(module, Acc1, Acc2, Histo,x1, x2, y1, y2, slope1, slope2);
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
    const char* filename = cimg_option("-i","/Users/rubcuevas/Desktop/Algorithmie de l'image/EllipseDetection/EllipseDetection/SharbatGula.bmp","Input image file");
    const double alpha   = cimg_option("-a",0.5,"Gradient smoothing");
    const double sigma   = cimg_option("-s",0.2,"Hough Transform smoothing");
    const int Rmin       = cimg_option("-rmin",5,"Minimum radius");
    const int Rmax       = cimg_option("-rmax",28,"Maximum radius");
    const int NbShape    = cimg_option("-n",20,"Number of the extracted parametric shapes");
    
    int ax;
 
 // Opening of filename
 CImg<> img(filename);
/*
 Added by me
 */
// Accumulator Computation
int NbTheta = 360;
int NbRho   = (int)(sqrt(img.width()*img.width()+img.height()*img.height())/2.);
CImg<> acc1(img.width(), img.height());
CImg<> acc2(NbTheta, NbTheta);
vector<int> Histo(img.width()/2);
    for(int i = 0; i < img.width()/2; ++i){
        Histo.push_back(0);
    }
// Display
CImgDisplay dispSpatial(img,"Input Image");
CImgDisplay acc1Spatial(acc1,"Accumulator map a0 and b0");
CImgDisplay acc2Spatial(acc2, "Accumulator KN");
 /* End of added by me */
 // Accumulator Computation
EllipseAccumulator(img.get_channel(0), acc1, acc2, Histo, alpha, sigma);
ax = FindMaxHisto(Histo);
 while (!dispSpatial.is_closed() && !acc1Spatial.is_closed())
 {
  dispSpatial.wait(dispSpatial,acc1Spatial);

  // When clicking on the vote window.
  if (acc1Spatial.button())
  {
   CImg<> win_param(1,3);
   win_param(0,0) = acc1Spatial.mouse_x();
   win_param(0,1) = acc1Spatial.mouse_y();
   DrawHoughLine(img,acc1,win_param,NbTheta,NbRho);

   dispSpatial.display(img);
   acc1Spatial.display(acc1);
  }
 }
 return 0;
}
