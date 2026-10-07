/* Axisymmetric Newtonian ideal MHD, first-order finite volume LLF + GLM.
   Units: L0, rho0, v0; mu0=1. Engine is prescribed, not a BH solution. */
class JetSolver {
 constructor(nr=56,nz=224,drive=1,ambient=1){this.nr=nr;this.nz=nz;this.h=12/nr;this.N=nr*nz;this.g=5/3;this.t=0;this.steps=0;this.drive=drive;this.ambient=ambient;this.on=true;this.floorE=0;this.injected=0;this.u=new Float64Array(this.N*9);this.du=new Float64Array(this.N*9);this.f=new Float64Array(9);this.fl=new Float64Array(9);this.fr=new Float64Array(9);this.ghost=new Float64Array(9);this.left=new Float64Array(9);this.right=new Float64Array(9);this.reset();}
 reset(){const {nr,nz,h,u}=this;for(let j=0;j<nz;j++)for(let i=0;i<nr;i++){let k=(j*nr+i)*9,r=(i+.5)*h,z=(j+.5-nz/2)*h;let rho=this.ambient,p=.08*this.ambient;u[k]=rho;u[k+6]=.18;u[k+4]=p/(this.g-1)+.5*.18**2;}this.initialEnergy=this.diagnostics().total;}
 pressure(a,k){return Math.max(1e-7,(this.g-1)*(a[k+4]-.5*(a[k+1]**2+a[k+2]**2+a[k+3]**2)/a[k]-.5*(a[k+5]**2+a[k+6]**2+a[k+7]**2)));}
 flux(a,k,d,out,ch){let rho=a[k],v=[a[k+1]/rho,a[k+2]/rho,a[k+3]/rho],b=[a[k+5],a[k+6],a[k+7]],p=this.pressure(a,k),bb=b[0]**2+b[1]**2+b[2]**2,pt=p+.5*bb,vn=v[d],bn=b[d],vb=v[0]*b[0]+v[1]*b[1]+v[2]*b[2];out[0]=rho*vn;for(let m=0;m<3;m++){out[1+m]=a[k+1+m]*vn-b[m]*bn+(m===d?pt:0);out[5+m]=m===d?a[k+8]:vn*b[m]-v[m]*bn;}out[4]=(a[k+4]+pt)*vn-vb*bn;out[8]=ch*ch*bn;let cs=this.g*p/rho,va=bb/rho;return Math.abs(vn)+Math.sqrt(.5*(cs+va+Math.sqrt(Math.max(0,(cs+va)**2-4*cs*bn*bn/rho))));}
 face(a,ka,b,kb,d,ch){
 const stride=d===0?9:this.nr*9,ia=ka/9%this.nr,ib=kb/9%this.nr;
 if(kb-ka===stride && ka-stride>=0 && kb+stride<b.length && (d===1||(ia>0&&ib<this.nr-1))){
 for(let m=0;m<9;m++){let dl=a[ka+m]-a[ka-stride+m],dr=a[ka+stride+m]-a[ka+m],sl=dl*dr>0?Math.sign(dl)*Math.min(Math.abs(dl),Math.abs(dr)):0;let el=b[kb+m]-b[kb-stride+m],er=b[kb+stride+m]-b[kb+m],sr=el*er>0?Math.sign(el)*Math.min(Math.abs(el),Math.abs(er)):0;this.left[m]=a[ka+m]+.5*sl;this.right[m]=b[kb+m]-.5*sr;}
 let valid=q=>q[0]>1e-5 && q[4]-.5*(q[1]**2+q[2]**2+q[3]**2)/q[0]-.5*(q[5]**2+q[6]**2+q[7]**2)>1e-7;
 if(valid(this.left)&&valid(this.right)){a=this.left;b=this.right;ka=0;kb=0;}
 }
 let s=Math.max(ch,this.flux(a,ka,d,this.fl,ch),this.flux(b,kb,d,this.fr,ch));for(let m=0;m<9;m++)this.f[m]=.5*(this.fl[m]+this.fr[m]-s*(b[kb+m]-a[ka+m]));}
 step(){let {u,du,nr,nz,h}=this;du.fill(0);let speed=.1;for(let k=0;k<u.length;k+=9){let rho=u[k],v=Math.hypot(u[k+1],u[k+2],u[k+3])/rho,b2=u[k+5]**2+u[k+6]**2+u[k+7]**2;speed=Math.max(speed,v+Math.sqrt((this.g*this.pressure(u,k)+b2)/rho));}let dt=.27*h/(2*speed),ch=speed;
 // Radial finite-volume areas / cylindrical volumes; axis flux area is zero.
 for(let j=0;j<nz;j++)for(let i=1;i<=nr;i++){let kl=(j*nr+i-1)*9,kr=i<nr?kl+9:kl;this.face(u,kl,u,kr,0,ch);for(let m=0;m<9;m++){du[kl+m]-=this.f[m]*i/((i-.5)*h);if(i<nr)du[kr+m]+=this.f[m]*i/((i+.5)*h);}}
 for(let j=0;j<=nz;j++)for(let i=0;i<nr;i++){let kl=(Math.max(0,j-1)*nr+i)*9,kr=(Math.min(nz-1,j)*nr+i)*9;this.face(u,kl,u,kr,1,ch);for(let m=0;m<9;m++){if(j>0)du[kl+m]-=this.f[m]/h;if(j<nz)du[kr+m]+=this.f[m]/h;}}
 for(let j=0;j<nz;j++)for(let i=0;i<nr;i++){let k=(j*nr+i)*9,r=(i+.5)*h,z=(j+.5-nz/2)*h,rho=u[k],vr=u[k+1]/rho,vp=u[k+3]/rho,br=u[k+5],bz=u[k+6],bp=u[k+7],pt=this.pressure(u,k)+.5*(br*br+bz*bz+bp*bp);
 du[k+1]+=(rho*vp*vp+pt-bp*bp)/r;du[k+3]+=(-rho*vr*vp+br*bp)/r;du[k+7]+=(vr*bp-vp*br)/r;du[k+5]+=u[k+8]/r;
 for(let m=0;m<9;m++)u[k+m]+=dt*du[k+m];u[k+8]*=Math.exp(-dt*ch/(.18*h));
 // Local electromotive forcing models winding by unresolved central engine.
 // Odd Bphi above/below plane; regular Bphi ~ r at axis; no axial momentum injected.
 let profile=r*Math.exp(-r*r/2.0)*Math.tanh(z/.45)*Math.exp(-z*z/1.5),ramp=1-Math.exp(-this.t/1.5),db=this.on?dt*8*this.drive*ramp*profile:0;
 let de=u[k+7]*db+.5*db*db;u[k+7]+=db;u[k+4]+=de;this.injected+=de*2*Math.PI*r*h*h;
 let dm=(this.on?.025:0)*dt*Math.exp(-r*r/2-z*z/1.5);u[k]+=dm;u[k+4]+=.08*dm/(this.g-1);
 if(u[k]<1e-5){u[k]=1e-5;}
 let min=.5*(u[k+1]**2+u[k+2]**2+u[k+3]**2)/u[k]+.5*(u[k+5]**2+u[k+6]**2+u[k+7]**2)+1e-7/(this.g-1);if(u[k+4]<min){this.floorE+=(min-u[k+4])*2*Math.PI*r*h*h;u[k+4]=min;}
 }this.t+=dt;this.steps++;return dt;}
 diagnostics(){let kin=0,mag=0,magT=0,thermal=0,maxV=0,div=0,total=0,head=0,flux=0;let {u,nr,nz,h}=this;for(let j=0;j<nz;j++)for(let i=0;i<nr;i++){let k=(j*nr+i)*9,r=(i+.5)*h,z=(j+.5-nz/2)*h,vol=2*Math.PI*r*h*h,rho=u[k],v2=(u[k+1]**2+u[k+2]**2+u[k+3]**2)/rho**2,b2=u[k+5]**2+u[k+6]**2+u[k+7]**2;kin+=.5*rho*v2*vol;mag+=.5*b2*vol;magT+=.5*u[k+7]**2*vol;thermal+=this.pressure(u,k)/(this.g-1)*vol;total+=u[k+4]*vol;maxV=Math.max(maxV,Math.sqrt(v2));if(Math.abs(u[k+7])>.12&&u[k+2]*z>0)head=Math.max(head,Math.abs(z));if(j===Math.floor(nz/2+4/h))flux+=u[k+2]*2*Math.PI*r*h; if(i>0&&i<nr-1&&j>0&&j<nz-1){let d=((r+h)*u[k+9+5]-(r-h)*u[k-9+5])/(2*h*r)+(u[k+nr*9+6]-u[k-nr*9+6])/(2*h);div=Math.max(div,Math.abs(d)*h/(Math.sqrt(b2)+1e-5));}}return {t:this.t,steps:this.steps,kin,mag,magT,thermal,total,maxV,div,head,flux,floorE:this.floorE, injected:this.injected};}
 snapshot(){let f=new Float32Array(this.N*8);for(let i=0;i<this.N;i++){let k=i*9,q=i*8,rho=this.u[k];f[q]=rho;f[q+1]=this.pressure(this.u,k);f[q+2]=this.u[k+1]/rho;f[q+3]=this.u[k+2]/rho;f[q+4]=this.u[k+5];f[q+5]=this.u[k+6];f[q+6]=this.u[k+7];f[q+7]=this.u[k+3]/rho;}return f;}
}
if(typeof module!=='undefined')module.exports=JetSolver;
if(typeof self!=='undefined'&&typeof document==='undefined'){let sim;self.onmessage=e=>{let m=e.data;if(m.type==='reset')sim=new JetSolver(m.nr||56,(m.nr||56)*4,m.drive,m.ambient);if(!sim)return;if(m.type==='step')for(let n=0;n<m.n;n++)sim.step();if(m.type==='engine')sim.on=m.on;let data=sim.snapshot();self.postMessage({data,diag:sim.diagnostics(),nr:sim.nr,nz:sim.nz,h:sim.h},[data.buffer]);};}
