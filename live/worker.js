let engine,probes=[];
self.onmessage=({data:m})=>{
 try{
  if(m.type==='init'){engine=new Circuit(m.data,m.params);probes=m.probes;if(m.demo)engine.pulse(m.stimulus,m.amplitude,m.duration,20);}
  if(m.type==='params')engine.configure(m.params);
  if(m.type==='probes')probes=m.probes;
  let pulse=null;if(m.type==='pulse')pulse=engine.pulse(m.index,m.amplitude,m.duration);
  const rows=[];if(m.type==='advance')for(let i=0;i<m.steps;i++){engine.step();rows.push([engine.time,...probes.flatMap(p=>engine.read(p))]);}
  const v=engine.v.slice(),pv=engine.pv.slice();
  self.postMessage({id:m.id,type:m.type,time:engine.time,params:engine.p,rows,v,pv,pulse},[v.buffer,pv.buffer]);
 }catch(e){self.postMessage({id:m.id,error:String(e.message||e)});}
};
