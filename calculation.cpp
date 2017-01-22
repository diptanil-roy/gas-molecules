#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <ctime>
#include <sstream>
#define pi 3.14159265

using namespace std;

time_t time_start= time(NULL);

const std::string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%Y_%m_%d.%X", &tstruct);
    return buf;
}

const int number_of_particles= 200;

ofstream fout("position.txt");
ofstream kout("numberofparticles.txt");
ofstream uout("speed.txt");
ofstream pout("pressure.txt");

class vector2d{
  public:
    double x,y, mag, angle;
    vector2d(){x=y=mag=angle=0;}
    vector2d(double a, double b){x=a; y=b; mag=abs(sqrt(pow(x,2) + pow(y,2))); angle=-atan2(y,x); }
    void fillvec(double a, double b){x=a; y=b; mag=sqrt(pow(x,2) + pow(y,2)); angle=-atan2(y,x); }
    vector2d operator + (vector2d v) {double xn, yn; xn=x+v.x; yn=y+v.y; vector2d vn(xn,yn); return vn;}
    vector2d operator - (vector2d v) {double xn, yn; xn=x-v.x; yn=y-v.y; vector2d vn(xn,yn); return vn;}
    vector2d operator - (){double xn,yn; xn=-x; yn=-y; vector2d vn(xn,yn); return vn;}
    vector2d operator * (double c) {double xn, yn; xn=x*c; yn=y*c; vector2d vn(xn,yn); return vn;}
    vector2d operator / (double);
    double operator * (vector2d v) {double dot= x*v.x+ y*v.y; return dot;}
    vector2d rotate(double theta) {double xn=x*cos(theta)+y*sin(theta); double yn=-x*sin(theta)+y*cos(theta); vector2d vn(xn,yn); return vn;}
};

vector2d vector2d::operator /(double c)
{
  double xn, yn;
   xn=x/c; yn=y/c; vector2d vn(xn,yn);
   return vn;
}

class particle{
  public:
    double x,y,ux,uy,size,mass; bool flag;
    particle(){x=y=ux=uy=size=0;}
    void fillval(double a, double b, double c, double d) {x=a; y=b; ux=c; uy=d; size=5; mass=1; flag=true;} //true: available for collision
    void newton2ndlaw(double, double);
    void boundary_check(double,double,double,double,double&, int&);
    friend int collision_check(particle, particle);
    vector2d* collision (particle, particle);
};

    double force_x(double x, double ux, double m)
    {
      double res= 0;
      return res;
    }

    double force_y(double y, double uy, double m)
    {
      double res=0;
      return res;
    }


    void particle::newton2ndlaw(double h, double m)
    {
      ux+= force_x(x,ux,m)*h;
      x+= ux*h;
      uy+= force_y(y,uy,m)*h;
      y+= uy*h;
    }

    void particle::boundary_check(double Dxt, double Dxb, double Dyt, double Dyb, double& impulse_per_length, int& boundary_collision_count)
    {
      if ((x+size)>=Dxt || (x+size)<= Dxb)
      {   ux *= -1; impulse_per_length+=2*mass*abs(ux)/(Dxt-Dxb); boundary_collision_count++;}
      else if ((y+size)>=Dyt || (y+size)<=Dyb)
      {   uy *= -1; }
    }

    int collision_check(particle p1, particle p2)
    {
      vector2d r1(p1.x,p1.y);
      vector2d r2(p2.x,p2.y);
      vector2d r=r1-r2;
      double distance_between_center = r.mag;
      if(distance_between_center <= p1.size+p2.size)
        return 1;
      else
        return 0;
    }

   vector2d* collision (particle a, particle b)
   {
     vector2d va2,vb2;
     vector2d zero_vector(0,0);
     vector2d va1(a.ux, a.uy); vector2d vb1(b.ux, b.uy);
     vector2d ra(a.x, a.y); vector2d rb(b.x,b.y);
     vector2d r=ra-rb;

      if(r.mag==0)
        {
          va2=va1;
          vb2=vb1;
        }
      else
        {
           double f = (va1-vb1)*r;
           va2= (va1 - r*(f/pow(r.mag,2)));       //Assuming collision is instantaneous.
           vb2= (vb1 + r*(f/pow(r.mag,2)));
        }

     vector2d* varr; varr= new vector2d[2];
     varr[0]=va2; varr[1]=vb2;

     double check_ec=va1.mag*va1.mag+vb1.mag*vb1.mag-va2.mag*va2.mag-vb2.mag*vb2.mag;
     vector2d check_pc=(va1+vb1)-(va2+vb2);

/*--------- energy and momentum conservation check------------------------------------
    if(abs(check_ec) < 0.5 && abs(check_pc.mag) < 0.5)
     {
       kout << "energy and momentum conserved \n";
     }
    else if (abs(check_ec) < 0.5)
    {
      kout << "energy conserved, but momentum not conserved\n";
    }
    else
    {
     kout << endl;
     kout << va1.x << "\t" << va1.y << "\t" << vb1.x << "\t" << vb1.y << endl;
     kout << va2.x << "\t" << va2.y << "\t" << vb2.x << "\t" << vb2.y << endl;
     kout << "Neither energy not conserved \n\n";
    }
//------------------------------------------------------------------------------------*/
     return varr;
   }

int main(){

      /*stringstream ss1, ss2, ss3, ss4;
      ss1 << "data/position_" << number_of_particles << "_" << currentDateTime() <<".txt" ;
      string s1 = ss1.str();
      const char* position_path=s1.c_str();
      fout.open(position_path);

      ss2 << "data/speed_" << number_of_particles << "_" << currentDateTime() <<".txt" ;
      string s2 = ss2.str();
      const char* speed_path=s2.c_str();
      uout.open(speed_path);

      ss3 << "debug/conservation_" << number_of_particles << "_" << currentDateTime() <<".txt" ;
      string s3 = ss3.str();
      const char* conservation_path=s3.c_str();
      kout.open(conservation_path);

      ss4 << "data/pressure_" << number_of_particles << "_" << currentDateTime() <<".txt" ;
      string s4 = ss4.str();
      const char* pressure_path=s4.c_str();
      pout.open(pressure_path);*/
      kout << number_of_particles;

      double x0, y0, ux0, uy0, t_final, Dxt, Dxb, Dyt, Dyb,h;
      t_final=120;
      Dxb = Dyb = 20; Dxt=800; Dyt=750;
      h=0.01;
      //srand (time(NULL));

      particle* parr= new particle[number_of_particles];

       for (int i=0; i<number_of_particles; i++)
         {
           x0=rand()%750+50;
           y0=rand()%650+50;
           ux0=rand()%500-250;
           uy0=rand()%500-250;
           parr[i].fillval(x0,y0,ux0,uy0);
        }

    double t, x_after_collision, y_after_collision, ux_after_collision, uy_after_collision, pressure, impulse_per_length, total_v2;
    int boundary_collision_count =0;

	cout << "Particles initialized. Calculation started. \n";

  for(t=0; t<=t_final; t=t+h)   //Master Loop: Drives time evolution through quantum h=1/60.
  {
    //---------Writing to File------------------------------------------------//
       fout << setprecision(3) << t;
       uout << setprecision(3) << t ;
       vector2d velocity;

      for(int i=0; i<number_of_particles; i++) //Write out the position of each particle in 2N columns.
        {
          fout << "\t" << setprecision(5) << parr[i].x << "\t" << setprecision(5) << parr[i].y << "\t";
          velocity.fillvec(parr[i].ux, parr[i].uy);
          uout << "\t" << setprecision(5) << velocity.mag  << "\t";   //For plotting M-B curve
        }
          uout << endl;
          fout << endl;
   //-------------------------------------------------------------------------//
      for(int i=0; i<number_of_particles; i++)
          {parr[i].flag=true; }                   //Make particles available for collision

      impulse_per_length=pressure=0;

      for(int i=0; i<number_of_particles; i++)          //Particle Loop: Place particle-specific routines here.
        {
          parr[i].newton2ndlaw(h,parr[i].mass);                      //Move the particle for next instant, just as required by the force.
          parr[i].boundary_check(Dxt, Dxb, Dyt, Dyb, impulse_per_length, boundary_collision_count);   //Check if it collides with the boundary.

          int j=i+1;
            while(j<number_of_particles && parr[i].flag)
             {
               if(parr[j].flag && collision_check(parr[i],parr[j])==1)
                       {
                         vector2d* varr=collision(parr[i],parr[j]);
                           x_after_collision=parr[i].x;
                           y_after_collision=parr[i].y;
                           ux_after_collision= varr[0].x;
                           uy_after_collision=varr[0].y;
                           parr[i].fillval(x_after_collision, y_after_collision, ux_after_collision, uy_after_collision);
                           parr[i].flag=false;

                           x_after_collision=parr[j].x;
                           y_after_collision=parr[j].y;
                           ux_after_collision= varr[1].x;
                           uy_after_collision=varr[1].y;
                           parr[j].fillval(x_after_collision, y_after_collision, ux_after_collision, uy_after_collision);
                           parr[j].flag=false;
                           j++;
                           break;
                       }
                  else
                      j++;
              }
          }               //Collisions for all the particles have been solved for this instant.

          total_v2=0;
          for(int i=0; i<number_of_particles; i++)
            {
              total_v2 += pow(parr[i].ux,2)+pow(parr[i].uy,2);
             }

          pressure=0.5*impulse_per_length/h;
          pout <<  setprecision(3) << t << "\t" << setprecision(5) << pressure  << endl;
         // kout << setprecision(4) << t << "\t" << setprecision(5) << total_v2  << endl;
     }  //Motion of the system has been solved for this instant.

  fout.close();
  uout.close();
  pout.close();


  time_t time_finish= time(NULL);
  //kout << "Total no of boundary collisions were:" << boundary_collision_count << '\n';
  kout.close();
  std::cout << "The program took " << time_finish-time_start << " seconds."<< '\n';
  return 0;
}
