
/*
My solution to 3d Convex Hull problem, 
For the internship application at WeRide.
This code is an implementation of the QuickHull algorithm (http://www.cogsci.rpi.edu/~destem/gamearch/quickhull.pdf).
By Puning Zhao
*/

#include<iostream>
#include<fstream>
#include<math.h>
#include<vector>
#include<string>

using namespace std;
int N;
vector< vector<double> > points;
const double eps=1e-8; //precision.

//Define 'vect' struct.
struct vect{
	double x,y,z;
};

/*
Now define 'face' struct. A face is determined by three vertices a,b,c. 
  The order of these three vertices are selected to ensure that crossproduct(findedge(a,b),findedge(a,c)) 
is always toward the outside of the polyhedron.
*/

struct face{
	int a,b,c;
	bool keep;
};

//Find the edge. The edge between point (x1,y1,z1) and point(x2,y2,z2) is defined as (x2-x1,y2-y1,z2-z1).
vect findedge(int i,int j){
	vect s;
	s.x=points[j][0]-points[i][0];
	s.y=points[j][1]-points[i][1];
	s.z=points[j][2]-points[i][2];
	return s;
}

//Find the square distance between two points.
double sqdist(int i,int j){
	return pow(points[j][0]-points[i][0],2)+pow(points[j][1]-points[i][1],2)+pow(points[j][2]-points[j][1],2);
}

//Define vector operations: dot product, cross product, vector sum.
double dotproduct(vect a,vect b){
	double result=a.x*b.x+a.y*b.y+a.z*b.z;
	return result;
}

vect crossproduct(vect a,vect b){
	vect v;
	v.x=a.y*b.z-a.z*b.y;
	v.y=a.z*b.x-a.x*b.z;
	v.z=a.x*b.y-a.y*b.x;
 return v;
}

vect vecsum(vect a, vect b){
	vect v;
	v.x=a.x+b.x;
	v.y=a.y+b.y;
	v.z=a.z+b.z;
	return v;
}

//vector length.
double length(vect v){
	return sqrt(dotproduct(v,v));
}

//Find area of the triangle constructed by three points.
double findarea(int i,int j,int k){
	vect edge1=findedge(i,j);
	vect edge2=findedge(i,k);
	double S=length(crossproduct(edge1,edge2))/2;
	return S;
}

//Find volume of the tetrahedron constructed by four points.
double findvolume(int i,int j, int k, int l){
	vect edge1=findedge(i,j);
	vect edge2=findedge(i,k);
	vect edge3=findedge(i,l);
	double V=dotproduct(edge3,crossproduct(edge1,edge2))/6;
	if(V<0) V=-V; //Find the absolute value of the volume.
	return V;
}

//Determine whether a face is visible to point i.
bool visible(int i,face f){
	vect edge1=findedge(f.a,f.b);
	vect edge2=findedge(f.a,f.c);
	double s=dotproduct(findedge(f.a,i),crossproduct(edge1,edge2));
	if (s>eps) return true;
	else return false;
}

class hull{
public:
	vector<face> faces;
	bool *conflict;
	int **indexf;
	int check(){
		int i,j;
		for(i=0;i<N;i++){
			if(conflict[i]){
				for(j=0;j<faces.size();j++){
					if(visible(i,faces[j])&&faces[j].keep) return i;
				}
			conflict[i]=false;
			} 
		}
		return -1;
	} //If there are still some points remaining, return the first index. If not, return -1.

	/*
	The following two functions (searchedge and searchreplace) performs deep first search to find all 
	the faces that are visible to point p. Note that when we incorporate p into our final convex hull,
	all the visible surfaces should be removed. It is also important to record the boundary of all 
	visible surfaces, so that we can add new surface by connecting point p with all the edges on the
	boundary.
	*/

	void searchedge(int p,int a, int b){
		int ind= indexf[a][b]; //Find the index of the new face.
		face newface;
		if(faces[ind].keep){
			if(visible(p,faces[ind])) {
				searchreplace(p,ind);
				if(p==8) cout<<"visible:"<<a<<b<<p<<endl;
			}
			else{
				/*
				If visible(p,faces[ind]) is false, then the edge (a,b) is an edge on the boundary of 
				all visible faces. The next few lines construct a new face connecting point p with 
				this edge.
				*/
				 newface.a=b;
				 newface.b=a;
				 newface.c=p;
				 newface.keep=true;
				 if(p==8) cout<<a<<" "<<b<<" "<<p<<endl;
				 faces.push_back(newface);
				 indexf[b][a]=faces.size()-1;
				 indexf[a][p]=faces.size()-1;
				 indexf[p][b]=faces.size()-1;
			}
		}
	}

	void searchreplace(int i,int j){ 
		if(i==8) cout<<j<<" "<<faces[j].a<<" "<<faces[j].b<<" "<<faces[j].c<<" "<<faces.size()<<endl;
		faces[j].keep=false;
		searchedge(i,faces[j].b,faces[j].a);
		searchedge(i,faces[j].c,faces[j].b);
		searchedge(i,faces[j].a,faces[j].c);
	}// This is a deep first search.

	void merge(int i){
		int M=faces.size();
		int j0,j,k;
		if(i==8){
			for(j=0;j<N;j++){
				for(k=0;k<N;k++){
					if(indexf[j][k]>=0) cout<<"  "<<indexf[j][k];
					else cout<<" "<<-1;
				}
				cout<<endl;
			}
		}
		for(j=0;j<M;j++){
			if(visible(i,faces[j])&&faces[j].keep){
				j0=j;
				break;
			}
		}
		cout<<"target: "<<j0<<endl;
	searchreplace(i,j0);
	conflict[i]=false;
	}

	void simplify(){
		vector<face> f;
		int i;
		for(i=0;i<faces.size();i++){
			if(faces[i].keep) f.push_back(faces[i]);
		}
		faces=f;
	}

	void output(){
		int i;
		for(i=0;i<faces.size();i++){
			cout<<faces[i].a<<" "<<faces[i].b<<" "<<faces[i].c<<endl;
		}
	}
};

hull initialize(){
	vector<face> Faces;
	hull Hull;
	int index[6]={0,0,0,0,0,0}; //Location of maximum x, minimum x, maximum y, minimum y, maximum z, minimum z.
	int p[4],i,j;
	double xmax,xmin,ymax,ymin,zmax,zmin,dist,S,V,maxdist=0,maxS=0,maxV=0;
	int t;
	face newface;

	Hull.conflict=new bool[N];
	Hull.indexf=new int* [N];
	for(i=0;i<N;i++){
		Hull.indexf[i]=new int[N];
	}

	/*
	The first step of QuickHull algorithm is to construct an initial hull using four points that do not fall on
	the same plane. The construction follows the following steps:
	1. Find the minimum and maximum value of each coordinate;
	2. Pick two points p1, p2 with maximum distance;
	3. Pick the third point, which maximizes the area of triangle p1-p2-p3;
	4. Pick the fourth point, which maximizes the volume of tetrahedron p1-p2-p3-p4;
	5. Construct an initial hull with p1-p2-p3-p4;
	6. Partite the points into different conflicting faces.
	*/

	xmax=xmin=points[0][0];
	ymax=ymin=points[0][1];
	zmax=zmin=points[0][2];
	for(i=1;i<N;i++){
		if(points[i][0]>xmax){
		xmax=points[i][0];
		index[0]=i;
		}
		if(points[i][0]<xmin){
			xmin=points[i][0];
			index[1]=i;
		}
		if(points[i][1]>ymax){
			ymax=points[i][1];
			index[2]=i;
		}
		if(points[i][1]<ymin){
			ymin=points[i][1];
			index[3]=i;
		}
		if(points[i][2]>zmax){
			zmax=points[i][2];
			index[4]=i;
		}
		if(points[i][2]<zmin){
			zmin=points[i][2];
			index[5]=i;
		}
	}
	p[0]=0;
	p[1]=0;
	for(i=0;i<6;i++){
		for(j=i+1;j<6;j++){
			dist=sqdist(index[i],index[j]);
			if(dist>maxdist){
				p[0]=index[i];
				p[1]=index[j];
				maxdist=dist;
			}
		}
	}
	for(i=0;i<N;i++){
		S=findarea(p[0],p[1],i);
		if(S>maxS){
			p[2]=i;
			maxS=S;
		}
	}
	for(i=0;i<N;i++){
		V=findvolume(p[0],p[1],p[2],i);
		if(V>maxV){
			p[3]=i;
			maxV=V;
		}
	}
	if(maxV<1e-6){
		cout<<"The points are almost on the same plane."<<endl;
		exit(-1);
	}
	cout<<"Initial Hull:"<<endl;
	cout<<p[0]<<" "<<p[1]<<" "<<p[2]<<" "<<p[3]<<endl;
	for(i=0;i<4;i++){
		newface.a=p[(i+1)%4];
		newface.b=p[(i+2)%4];
		newface.c=p[(i+3)%4];
		if(visible(p[i],newface)){
			t=newface.b;
			newface.b=newface.c;
			newface.c=t;
		} 
		/*
		Two comments:
		1.Recall that we need to ensure that the crossproduct(findedge(a,b),findedge(a,c)) 
		is always toward the outside of the polyhedron. Here I swap b,c if this condition
		is not satisfied, i.e. newface is visible from point p[i].
		
		2.The QuickHull algorithm requires us to find the neighbor of the surface. Therefore,
		I created a two dimensional array indexf, such that indexf[a][b] and indexf[b][a] 
		represent two faces which intersects at line (a,b) separately.
		*/
		Faces.push_back(newface);
		Hull.indexf[newface.a][newface.b]=Faces.size()-1;
		Hull.indexf[newface.b][newface.c]=Faces.size()-1;
		Hull.indexf[newface.c][newface.a]=Faces.size()-1;
	}
	Hull.faces=Faces;
	for(i=0;i<=N;i++){
		Hull.conflict[i]=true;
	}
	for(i=0;i<4;i++){
		Hull.conflict[p[i]]=false;
	}
	return Hull;
}

//main function
int main(){
	int i;
	hull Hull;
	double a,b,c;
	char filename[50];
	ifstream file;
	cout<<"Enter file name:"<<endl;
	cin>>filename;
	file.open(filename);
	if(!file){
		cout<<"File not found"<<endl;
		exit(-1);
	}
	while(file>>a){
		vector<double> temp;
		temp.push_back(a);
		file>>b;
		file>>c;
		temp.push_back(b);
		temp.push_back(c);
		points.push_back(temp);
	}
	N=points.size();
	cout<<"Number of points: "<<N<<endl;
	if(N<4){
		cout<<"The number of points should be at least 4."<<endl;
		exit(-1);
	}
	Hull=initialize();
	i=Hull.check();
	while(i!=-1){
		cout<<"Merging point "<<i<<endl;
		Hull.merge(i);
		cout<<"Merging complete"<<endl;
		i=Hull.check();
	} 
	Hull.simplify();
	cout<<"Result:"<<endl;
	Hull.output();
	return 0;
}