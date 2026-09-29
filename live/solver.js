/* GPL-3.0. Full-resolution forest elimination; no morphology downsampling. */
const DEFAULTS={radius:.25,resistivity:150,cm:1,tau:20,cellC:.5,gain:.05,scale:5,sign:1,axial:1,dt:.5,feedback:true,cut:-1};
class Circuit {
 constructor(D,p={}) {
  this.D=D;this.n=D.baseC.length;this.time=0;this.steps=0;this.v=new Float64Array(this.n);this.pv=new Float64Array(D.cells.length);
  this.parent=new Int32Array(this.n).fill(-1);this.parentEdge=new Int32Array(this.n).fill(-1);this.order=[];
  const adj=Array.from({length:this.n},()=>[]);D.edges.forEach(([a,b],k)=>{adj[a].push([b,k]);adj[b].push([a,k])});
  const seen=new Uint8Array(this.n);
  for(let root=0;root<this.n;root++)if(!seen[root]){seen[root]=1;const q=[root];for(let at=0;at<q.length;at++){const i=q[at];this.order.push(i);for(const [j,k] of adj[i])if(!seen[j]){seen[j]=1;this.parent[j]=i;this.parentEdge[j]=k;q.push(j)}else if(j!==this.parent[i])throw Error('Axial anatomy must be a forest');}}
  this.rhs=new Float64Array(this.n);this.prhs=new Float64Array(D.cells.length);this.pulses=[];this.configure({...DEFAULTS,...p});
 }
 configure(p) {
  const next={...(this.p||DEFAULTS),...p};
  for(const k of ['radius','resistivity','cm','tau','cellC','scale','dt'])if(!Number.isFinite(next[k])||next[k]<=0)throw Error(k+' must be positive');
  for(const k of ['gain','axial'])if(!Number.isFinite(next[k])||next[k]<0)throw Error(k+' must be nonnegative');
  if(![-1,1].includes(next.sign)||!Number.isInteger(next.cut)||next.cut < -1||next.cut>=this.D.edges.length)throw Error('Invalid sign or cut');
  if(this.pulses.some(x=>Math.abs(x.start/next.dt-Math.round(x.start/next.dt))>1e-7||Math.abs(x.end/next.dt-Math.round(x.end/next.dt))>1e-7))throw Error('Pulse boundaries must align with dt');
  this.p=next;const D=this.D;
  this.C=Float64Array.from(D.baseC,(c,i)=>D.xyz[i]?c*(next.radius/.25)*next.cm:next.cellC);
  this.pc=new Float64Array(D.cells.length);this.C.forEach((c,i)=>this.pc[D.owner[i]]+=c);
  this.diag=Float64Array.from(this.C,c=>c/next.dt+c/next.tau);this.g=new Float64Array(this.n);
  for(const i of this.order){const par=this.parent[i],e=this.parentEdge[i];if(par>=0){const g=e===next.cut?0:D.edges[e][2]*(next.radius/.25)**2*(150/next.resistivity)*next.axial;this.g[i]=g;this.diag[i]+=g;this.diag[par]+=g;}}
  for(let k=this.order.length-1;k>=0;k--){const i=this.order[k],par=this.parent[i];if(par>=0)this.diag[par]-=this.g[i]**2/this.diag[i];}
  for(const x of this.diag)if(!Number.isFinite(x)||x<=0)throw Error('Invalid passive solve');
 }
 pulse(index,amplitude,duration,start=this.time){
  if(!Number.isInteger(index)||index<0||index>=this.n||!Number.isFinite(amplitude)||!Number.isFinite(duration)||duration<=0||!Number.isFinite(start)||start<this.time-1e-7)throw Error('Invalid pulse');
  const dt=this.p.dt;start=Math.round(start/dt)*dt;const end=start+Math.max(1,Math.round(duration/dt))*dt;
  this.pulses.push({index,amplitude,start,end});return this.pulses[this.pulses.length-1];
 }
 step(){
  const {D,p,v,pv,rhs,prhs}=this;
  for(let i=0;i<this.n;i++)rhs[i]=this.C[i]/p.dt*v[i];
  for(let c=0;c<pv.length;c++)prhs[c]=this.pc[c]/p.dt*pv[c];
  for(let k=0;k<D.pre.length;k++)if(p.feedback||!D.feedback[k]){const a=D.pre[k],b=D.post[k];rhs[b]+=p.sign*p.gain*Math.tanh(v[a]/p.scale);prhs[D.owner[b]]+=p.sign*p.gain*Math.tanh(pv[D.owner[a]]/p.scale);}
  for(const x of this.pulses)if(this.time>=x.start-1e-8&&this.time<x.end-1e-8){rhs[x.index]+=x.amplitude;prhs[D.owner[x.index]]+=x.amplitude;}
  for(let k=this.order.length-1;k>=0;k--){const i=this.order[k],j=this.parent[i];if(j>=0)rhs[j]+=this.g[i]/this.diag[i]*rhs[i];}
  for(const i of this.order){const j=this.parent[i];v[i]=(rhs[i]+(j>=0?this.g[i]*v[j]:0))/this.diag[i];if(!Number.isFinite(v[i]))throw Error('Nonfinite voltage; reset with smaller gain or timestep');}
  for(let c=0;c<pv.length;c++)pv[c]=prhs[c]/(this.pc[c]/p.dt+this.pc[c]/p.tau);
  this.time=Number((this.time+p.dt).toFixed(9));this.steps++;
  this.pulses=this.pulses.filter(x=>x.end>this.time+1e-8);
 }
 read(probe){
  const {D}=this;if(probe.kind==='node')return [this.v[probe.index],this.pv[D.owner[probe.index]]];
  const indices=probe.kind==='cell'?D.cells[probe.index].indices:D.cells.filter(c=>c.type===probe.type).flatMap(c=>c.indices);
  let s=0,point=0,cap=0;for(const i of indices){s+=this.C[i]*this.v[i];point+=this.C[i]*this.pv[D.owner[i]];cap+=this.C[i];}return [s/cap,point/cap];
 }
}
if(typeof module!=='undefined')module.exports={Circuit,DEFAULTS};
