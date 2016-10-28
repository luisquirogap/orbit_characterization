
#include<../useful_routines/useful_routines.c>

int main(int argc, char *argv[])
{
  int i, nPoints, counter;
  int indexMin, indexMax, nPartDisk;
  int nMin, nMax;
  double tf, rf, vf;
  char *infile, namefiles[200];

  CM cmDisk;
  
  FILE *fOrbit, *fInter, *fAngMom, *fPeriastro;
  
  infile = argv[1];
  indexMin = atoi(argv[2]);
  indexMax = atoi(argv[3]);
  nPartDisk = atoi(argv[4]);
  
  sprintf(namefiles,"%s_satOrbit.output",infile);
  
  fOrbit = fopen(namefiles,"r");
  
  nPoints = counterLines(namefiles);
  
  int time[nPoints];
  double x[nPoints], y[nPoints], z[nPoints];
  double vx[nPoints], vy[nPoints], vz[nPoints];
  double times[nPoints], r[nPoints], v[nPoints], tq, q, vq;
  double tData, qData, vData, L[3], LMag, LDisk;
  double LDot, inclination, meanInclination, iniInclination, periInclination;
  
  sprintf(namefiles,"%s_orbital_momentum.output",infile);
  fAngMom = fopen(namefiles,"w");
  
  qData = 1e8;
  meanInclination = 0.0;
  counter = 0;
  
  for( i=0; i<nPoints; i++)
    {
  
      returnRead = fscanf(fOrbit,"%d %lf %lf %lf %lf %lf %lf",&time[i],&x[i],&y[i],&z[i],&vx[i],&vy[i],&vz[i]);
      times[i] = 1.0*time[i];
      r[i] = sqrt( x[i]*x[i] + y[i]*y[i] + z[i]*z[i] );
      v[i] = sqrt( vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i] );
      
      L[X] = y[i]*vz[i] - z[i]*vy[i];
      L[Y] = z[i]*vx[i] - x[i]*vz[i];
      L[Z] = x[i]*vy[i] - y[i]*vx[i];
      
      LMag = sqrt( L[X]*L[X] + L[Y]*L[Y] + L[Z]*L[Z] );
      
      sprintf(namefiles,"%s_%.3d",infile,i);
      
      read_gadget1(namefiles);
      
      nMin = N_part[0] + N_part[1];
      nMax = nMin + nPartDisk;
      
      compute_angmom(NULL, nPartDisk, nMin, nMax, &cmDisk);
      printf("LCM : %e %e %e\n",cmDisk.lcm[X],cmDisk.lcm[Y],cmDisk.lcm[Z]);
      
      LDisk = sqrt( cmDisk.lcm[X]*cmDisk.lcm[X] +
		    cmDisk.lcm[Y]*cmDisk.lcm[Y] +
		    cmDisk.lcm[Z]*cmDisk.lcm[Z] );
      
      LDot = L[X]*cmDisk.lcm[X] + L[Y]*cmDisk.lcm[Y] + L[Z]*cmDisk.lcm[Z];
      
      inclination = LDot/(LMag*LDisk);
      
      inclination =  acos(inclination)*180.0/M_PI;
      
      if( i==0 )
	iniInclination = inclination;
      
      fprintf(fAngMom,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
	      times[i],0.0
	      ,L[X],L[Y],L[Z],LMag,
	      cmDisk.lcm[X],cmDisk.lcm[Y],cmDisk.lcm[Z],LDisk,
	      inclination);

      if( ( i>=indexMin ) && ( i<=indexMax) )
	if( r[i] < qData )
	  {
	    qData = r[i];
	    vData = v[i];
	    tData = times[i];
	    periInclination = inclination;
	   
	  }

      meanInclination = meanInclination + inclination;
      counter++;
      
      free(particles);

    }

  //  meanInclination = meanInclination/(indexMax-indexMin);
  // meanInclination = meanInclination/nPoints;
  meanInclination = meanInclination/counter;
  
  fclose(fOrbit);
  fclose(fAngMom);
  
  sprintf(namefiles,"%s_interpolated_orbit.output",infile);
  fInter = fopen(namefiles,"w");
  
  gsl_interp_accel *accr = gsl_interp_accel_alloc();
  gsl_interp_accel *accv = gsl_interp_accel_alloc();
  gsl_spline *spliner = gsl_spline_alloc(gsl_interp_cspline, nPoints);
  gsl_spline *splinev = gsl_spline_alloc(gsl_interp_cspline, nPoints);
  
  gsl_spline_init (spliner, times, r, nPoints);
  gsl_spline_init (splinev, times, v, nPoints);
  
  q = r[0];
  
  for( tf = times[indexMin]; tf<times[indexMax]; tf += 0.001 )
    {
      rf = gsl_spline_eval(spliner, tf, accr);
      vf = gsl_spline_eval(splinev, tf, accv);
      
      if( rf < q )
	{
	  tq = tf;
	  q = rf;
	  vq = vf;
	}
      
      fprintf(fInter,"%lf %lf %lf\n",tf,rf,vf);
    }
  
  sprintf(namefiles,"%s_periastro.output",infile);
  fPeriastro = fopen(namefiles,"w");
  fprintf(fPeriastro,"%lf %lf %lf\n", tq, q, vq);
  fprintf(fPeriastro,"%lf %lf %lf\n", tData, qData, vData);
  fprintf(fPeriastro,"%lf %lf %lf\n", iniInclination, periInclination, meanInclination);
  
  gsl_spline_free (spliner);
  gsl_interp_accel_free (accr);
  gsl_spline_free (splinev);
  gsl_interp_accel_free (accv);
  
  fclose(fInter);
  fclose(fPeriastro);
  
  return 0;
}

