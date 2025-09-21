//**********************************
// StochasticFractal.java
// author: Non-Euclidean Dreamer
//**********************************

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Random;

import javax.imageio.IIOImage;
import javax.imageio.ImageIO;

public class StochasticFractal 
{
	private static String name="M",
			affine="affine";
	private static int width=2560,height=1440,depth=500, it=500000,start=9600;
	static int[]dimensions= {width,height,depth};
	static BufferedImage image=new BufferedImage(width,height,BufferedImage.TYPE_4BYTE_ABGR);
	static int dim=2;
	int cl=512*325;
	static int color=Color.green.getRGB();
	static int[] c= {0,255,255},dir= {1,-1,1};
	static double third=1.0/3,zoom;
	static double[]focus= {-.03,0.47},min={0,0},max={0,0};
	String type=affine;
	Random rand;
	static DecimalFormat df=new DecimalFormat("0000");
	double[][]parameters;
	double[]probabilities;
	 
	public static void main(String[]args)
	{
		int[]jc= {1,-1,0,0,0,0,0,0,0,0,-1,1,1,-1,-1,-1,1,1,1,1,1,1,-1,-1,-1,-1},
			 kc= {0,0,-1,1,0,0,-1,1,1,-1,0,0,0,0,-1,1,1,-1,-1,-1,1,1,1,1,-1,-1},
			 lc= {0,0,0,0,1,-1,-1,-1,1,1,1,1,-1,-1,0,0,0,0,1,-1,-1,1,1,-1,-1,1};
		zoom=704;//509.2967137900334{0.4944017692309294,0.4888529618459649,}M2700.png finished. 605.6242471319865{0.4990737286749461,0.49818429524319174,}M2880.png 528.0859625631399{0.5000580405614472,0.4951134838011577,}M3060.png
		StochasticFractal fractal=random(4);//new StochasticFractal(affine,new double[][]{ {cos0*.5,.5*sin0,-.5+cos0*.5,-.5*sin0,0.5*cos0,1+sin0*.5  }, {cos1*.5,0.5*sin1,0,-.5*sin1,0.5*cos1,1 }, {cos2*.5,0.5*sin2,.5*(1-cos2),-.5*sin2,0.5*cos2,1-0.5*sin2 }, },//{{.5*cos,-.5*sin,1,.5*sin,.5*cos,0},{.5,0,.5,0,.5,sqrt*.5},{.5*cos,.5*sin,1.5,-.5*sin,.5*cos,sqrt*.5}},//{-0.015252008272264028, 0.10372902577437082, 0.28143872542791204, 0.9356596848665755, 0.5229587681514538, 1.2936713932049626, }, 	{-0.8810553836675441, -0.5132364987897058, 0.6411026689004511, 0.9108058289663559, -0.44696471106126845, -0.0067947273036477185, } },
		double[] v= {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
		loadImage("M"+df.format(start-1)+".png");//
int t=12;for(int i=0;i<720;i++) //	int i=360;for(int t=0;t<jc.length;t++)//				
		{int j=-360,k=720,l=360,m=720;
		double q=1-i/360.0,r=-1;
		
		double sqrt=Math.sqrt(3),cos0=Math.cos(Math.PI*(j/720.0)),sin0=Math.sin(Math.PI*(j/720.0)),cos1=Math.cos(Math.PI*(k/720.0)),sin1=Math.sin(Math.PI*(k/720.0)),cos2=Math.cos(Math.PI*(l/720.0)),sin2=Math.sin(Math.PI*(l/720.0)),cos3=Math.cos(Math.PI*(m/720.0)),sin3=Math.sin(Math.PI*(m/720.0));
			//new double[]{.3,0.3,0.4}) ;
		fractal=new StochasticFractal(affine,new double[][]{ {cos0*1,sin0*0,-1,-sin0*q,cos0*-0,.5}, {cos1*(1.0/3),sin1*(1.0/6),.5,-sin1*2.0/3*q,0.5*cos1*q,1}, {cos2*.5,sin2/3,1.0/3,-.5*sin2*q,cos2/3*q,.5},{cos3*-(1.0/6),sin3*0,0,-sin3*.5*q,cos3/-3*q,.5*(1-sin2)} },//{{.5*cos,-.5*sin,1,.5*sin,.5*cos,0},{.5,0,.5,0,.5,sqrt*.5},{.5*cos,.5*sin,1.5,-.5*sin,.5*cos,sqrt*.5}},//{-0.015252008272264028, 0.10372902577437082, 0.28143872542791204, 0.9356596848665755, 0.5229587681514538, 1.2936713932049626, }, 	{-0.8810553836675441, -0.5132364987897058, 0.6411026689004511, 0.9108058289663559, -0.44696471106126845, -0.0067947273036477185, } },
				new double[] {.25,.25,.25,.25});

		dampen(0.9);//
		//fractal.print();
		//color=spectrum((i+start)*10,1,1);
		fractal.draw(new double[] {0,0},  it);
	//	fractal.parameters[4][0]+=1.0/128
		print(name+df.format(start+i)+".png");
		
		
	//	c[0]++;
	//	color=new Color(c[0],c[1],c[2]).getRGB();
	//	fractal.edit(v);print(v);
	//	image=new BufferedImage(1080,1080,BufferedImage.TYPE_4BYTE_ABGR);
	}
		
	}
	
	private void shift(int i,double[]v)
	{
			for(int j=0;j<dim;j++)
			parameters[i][j*(dim+1)+dim]+=v[j];
	}
	
	
	
	private static void print(double[] v) {
		System.out.print("{");
		for(int i=0;i<v.length;i++)
			System.out.print(v[i]+",");
		System.out.print("}");
	}


	private static void dampen(double d) {
		for(int i=0;i<image.getWidth();i++)
			for(int j=0;j<image.getHeight();j++)
				image.setRGB(i, j,  darken(image.getRGB(i, j),d));
	}


	private static int darken(int rgb, double d) {
		Color c=new Color(rgb);
		int r=(int) (c.getRed()*d), g=(int) (c.getGreen()*d),b=(int) (c.getBlue()*d);
		return new Color(r,g,b).getRGB();
	}


	private void edit(double[] v) {
		int k=rand.nextInt(v.length);
		
		v[k]+=.001*(rand.nextDouble()-.5);
		
		for(k=0;k<v.length;k++)v[k]*=0.9999;
		double bound;
		double norm,sum=0;
		double[]det=new double[parameters.length];
		
		for(int i=0;i<parameters.length;i++)
		{	for(int j=0;j<parameters[0].length;j++)
			{
				parameters[i][j]+=v[i*parameters.length+j];
				norm=Math.abs(parameters[i][j]);
				bound=.99;
				if(j%(dim+1)==dim) {bound=2;if(norm>bound)parameters[i][j]*=bound/norm;}
				if(j%(dim+1)==dim-1) {for(int l=1;l<dim;l++)norm+=Math.abs(parameters[i][j-l]);
							if(norm>bound) {for(int l=0;l<dim;l++)parameters[i][j-l]*=bound/norm;}}
				
			
			}
			det[i]=Math.abs(det(parameters[i]));
			sum+=det[i];
		}
		for(int i=0;i<probabilities.length;i++)probabilities[i]=det[i]/sum;
	}
	
	private static double det(double[] ds) {
		if(dim==2)
		return ds[0]*ds[4]-ds[1]*ds[3];
		
		return ds[0]*ds[5]*ds[10]+ds[1]*ds[6]*ds[8]+ds[2]*ds[4]*ds[9]-ds[2]*ds[5]*ds[8]-ds[0]*ds[6]*ds[9]-ds[1]*ds[4]*ds[10];
	}


	public static void loadImage(String name)
	{
		File file=new File(name);
		
		try {
			image=ImageIO.read(file);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}


	private void print() {
		System.out.print("{");
		for(int i=0;i<probabilities.length;i++)System.out.print(probabilities[i]+", ");
		System.out.print("} - { ");
		for(int i=0;i<probabilities.length;i++) {System.out.print("{");
		for(int j=0;j<6;j++)
			System.out.print(parameters[i][j]+", ");
		System.out.print("}, ");
		}
		System.out.println("}");
	}


	private static void print(String n) {
		File outputfile = new File(n);
		try {  
			ImageIO.write(image,"png", outputfile);
			System.out.println(n+" finished.");
		} catch (IOException e) {
			System.out.println("IOException");
			e.printStackTrace();
		}
	}


	public StochasticFractal(String typ, double[][]par, double[]p)
	{
		type=typ;
		parameters=par;
		probabilities=p;
		rand=new Random();
	}
	
	public void draw(double[] start, int it)
	{
		double[] point=start.clone();
		cl=0;
		for(int i=0;i<it;i++)
		{
			draw(point,focus,zoom);
			iterate(point);
	//		print(point);	System.out.println();
		}
		
	//	print(min);print(max);
		zoom=Math.max (1,Math.pow(zoom,.99)*Math.pow(.95*Math.min(width/(max[0]-min[0]), height/(max[1]-min[1])),.01));
		if(dim>2)zoom=Math.max (1,Math.pow(zoom,.99)*Math.pow(.95*Math.min(width/(max[0]-min[0]), height/(max[1]-min[1])),.01));
		for(int i=0;i<2;i++)
		{
			focus[i]=focus[i]*.99+.005*(max[i]+min[i]);
			max[i]=0;min[i]=0;
		}
		if(dim>2) {focus[2]=focus[2]*.99+.005*(3*min[2]-max[2]);max[2]=0;min[2]=0;}
		System.out.println();
		System.out.print(zoom);
		print(focus);
	}

	private void iterate(double[] point) 
	{
		double[]nju=new double[dim];
		double r=rand.nextDouble();
		int k=0,d1=dim+1;
		while(r>probabilities[k]) 
		{
			r-=probabilities[k];
			k++;
		}
		cl=(cl+ 6*256*k)/probabilities.length;color=spectrum(cl,	1,1);

		if(type==affine) {
			for(int i=0;i<dim;i++)
			{
				for(int j=0;j<dim;j++)
				{
					nju[i]+=point[j]*parameters[k][i*d1+j];
				}
				nju[i]+=parameters[k][i*d1+dim];
			}
		}
		for(int i=0;i<dim;i++) {point[i]=nju[i];	min[i]=Math.min(min[i], point[i]);max[i]=Math.max(max[i], point[i]);}
	
	}

	private static void draw(double[] point, double[] loc, double zoom) {
		int x=(int) ((point[0]-loc[0])*zoom+width/2),
				y=(int) ((-point[1]+loc[1])*zoom+height/2);
	
		try {image.setRGB(x, y, color);
	}catch(IndexOutOfBoundsException e) {}
}
	private static StochasticFractal random(int n)
	{
		Random rand=new Random();
		double[] prob=new double[n];
		double[][]par=new double[n][6];
		for(int i=0;i<n;i++)
		{
			prob[i]=1.0/n;
			for(int j=0;j<6;j++)
			{
				par[i][j]=(rand.nextDouble()-.5)*2;
				if(j%3==1)par[i][j]*=1-Math.abs(par[i][j-1]);
				if(j%3==2)par[i][j]*=2;
			}
		}
		return new StochasticFractal(affine,par,prob);
	}
	
	
	private static StochasticFractal Barnsley() {
		double[]prob= {0.01,0.85,0.07,0.07};
		double[][]par= {{0,0,0,0,0.16,0},{0.85,0.04,0,-.04,.85,1.6},{.2,-.26,0,.23,.22,1.6},{-.15,.28,0,.26,.24,.44}
		};
		return new StochasticFractal(affine,par,prob);
	}
	private static StochasticFractal Koch(double d) {
		double[]prob= {0.5,0.5};
		double phi=Math.PI*d/180, l=-.5/Math.cos(phi), cos=Math.cos(phi)*l,sin=Math.sin(phi)*l;
		double[][]par= {{-cos,sin,0,sin,cos,0},{cos,-sin,1,sin,cos,0}};
		
		return new StochasticFractal(affine,par,prob);
	}
	public static int spectrum(int n, double d, double l)
	{
		double full=d*l, term=(1-l)*(1+d)*255/2;
		n=(n%(6*256)+6*256)%(6*256);
		if (n<256)
			return new Color((int) (255*full+term),(int) (n*full+term),(int) term).getRGB();
		n-=256;
		if (n<256)
			return new Color((int) ((255-n)*full+term),(int) (255*full+term),(int) term).getRGB();
		n-=256;
		if (n<256)
			return new Color((int) term,(int) (255*full+term),(int) (n*full+term)).getRGB();
		n-=256;
		if (n<256)
			return new Color((int) term,(int) ((255-n)*full+term),(int) (255*full+term)).getRGB();
		n-=256;
		if (n<256)
			return new Color((int) (n*full+term),(int) term,(int) (255*full+term)).getRGB();
		n-=256;
			return new Color((int) (255*full+term),(int) term,(int) ((255-n)*full+term)).getRGB();
	}
	
	static double[][][]pars={};
}
