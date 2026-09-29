const fs=require('fs');const {Circuit}=require('../live/solver.js');const D=JSON.parse(fs.readFileSync(process.argv[2]));const cases=JSON.parse(fs.readFileSync(process.argv[3]));let out=[];
for(const c of cases){const e=new Circuit(D,c.params);e.pulse(c.stimulus,c.amplitude,c.duration,20);for(let i=0;i<320;i++){if(c.change&&i===80)e.configure(c.change);e.step();}out.push({name:c.name,spatial:Array.from(e.v),point:Array.from(e.pv),time:e.time});}
fs.writeFileSync(process.argv[4],JSON.stringify(out));
