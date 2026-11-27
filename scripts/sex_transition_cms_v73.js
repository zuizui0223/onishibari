// Hypothetical cytoplasmic male sterility and autosomal dominant pistil-loss modifier.
function run(){
 const mend=(a,b)=>[(1-a)*(1-b),a*(1-b)+(1-a)*b,a*b];
 const female=i=>i===3;
 const seeds=(i,p)=>female(i)?p.A:i>=3?p.w:1;
 const pollen=(i,p)=>female(i)?0:(i%3!==0&&(p.mode==="constitutive"||i<3)?1-p.k:1);
 const S=x=>x.reduce((a,b)=>a+b,0);
 function resident(x,p){
   const W=S(x.map((v,i)=>v*pollen(i,p))); if(W<1e-13)return null;
   const r=S(x.map((v,i)=>v*pollen(i,p)*(i%3)/2))/W,T=S(x.map((v,i)=>v*seeds(i,p)));
   const y=Array(6).fill(0);
   for(let i=0;i<6;i++){
     const c=3*Math.floor(i/3),a=(i%3)/2,t=female(i)?0:p.t;
     const o=mend(a,r),self=mend(a,a);
     for(let j=0;j<3;j++)y[c+j]+=x[i]*seeds(i,p)*((1-t)*o[j]+t*self[j])/T;
   }
   if(Math.abs(S(y)-1)>1e-10||y.some(v=>v<0))throw Error("resident mass");
   return y;
 }
 function rare(x,z,p,g){
   const W=S(x.map((v,i)=>v*pollen(i,p))),T=S(x.map((v,i)=>v*seeds(i,p)));
   const r=S(x.map((v,i)=>v*pollen(i,p)*(i%3)/2))/W;
   // Rare heterozygous modifier is expressed in H only, zero maternal recruitment there.
   const y=Array(6).fill(0);
   const f=mend(0,r);
   for(let j=0;j<3;j++)y[3+j]+=z[3]*p.A*.5*f[j]/T;
   const rareP=S(z.map((v,i)=>v*pollen(i,p)*(1+g)))/W;
   const rareR=S(z.map((v,i)=>v*pollen(i,p)*(1+g)*(i%3)/2))/W;
   for(let i=0;i<6;i++){
     const c=3*Math.floor(i/3),a=(i%3)/2,o=x[i]*seeds(i,p)*(female(i)?1:1-p.t);
     for(let j=0;j<3;j++){
       const r0=mend(a,0)[j],r1=mend(a,1)[j];
       y[c+j]+=.5*o*((rareP-rareR)*r0+rareR*r1)/T;
     }
   }
   return y;
 }
 function full(x,p,g){
   // Explicit cytotype x restorer genotype x modifier genotype. Used for finite-difference test.
   const w=i=>pollen(Math.floor(i/3),p)*(i%3===0?1:1+g);
   const F=i=>female(Math.floor(i/3));
   const m=i=>F(i)?p.A:i%3===0?seeds(Math.floor(i/3),p):0;
   const W=S(x.map((v,i)=>v*w(i))),T=S(x.map((v,i)=>v*m(i)));
   const gam=[0,0,0,0];
   function gp(i){let a=(Math.floor(i/3)%3)/2,b=(i%3)/2;return [(1-a)*(1-b),(1-a)*b,a*(1-b),a*b]}
   for(let i=0;i<18;i++){let h=gp(i);for(let j=0;j<4;j++)gam[j]+=x[i]*w(i)*h[j]/W}
   const y=Array(18).fill(0);
   for(let i=0;i<18;i++){let c=9*Math.floor(i/9),t=F(i)?0:p.t,h=gp(i);
     for(let u=0;u<4;u++)for(let v=0;v<4;v++){
       let r=(u>>1)+(v>>1),b=(u%2)+(v%2);
       y[c+3*r+b]+=x[i]*m(i)*h[u]*((1-t)*gam[v]+t*h[v])/T;
     }
   }return y;
 }
 let checks=0;
 for(let mode of ["constitutive","silent"])for(let t of [0,.4])for(let g of [0,.7]){
   let p={A:1.5,w:.7,k:.3,mode,t},x=[.23,.17,.1,.2,.2,.1],z=[.1,.2,.1,.2,.3,.1],eps=1e-8;
   let v=Array(18).fill(0);for(let i=0;i<6;i++){v[3*i]=(1-eps)*x[i];v[3*i+1]=eps*z[i]}
   let y=full(v,p,g),r=rare(x,z,p,g),a=resident(x,p);
   if(Math.abs(S(y)-1)>1e-10||r.some((v,i)=>Math.abs(y[3*i+1]/eps-v)>1e-7)||a.some((v,i)=>Math.abs(y[3*i]-v)>1e-7))throw Error("full recurrence check");checks++;
 }
 function trajectory(p,start){
   let x=start,seq=[],n=4000;
   for(let i=0;i<n;i++){x=resident(x,p);if(!x)return {status:"pollen_collapse"};if(i>=n-300)seq.push(x)}
   let qs=seq.map(v=>v[3]),cs=seq.map(v=>S(v.slice(3))),q=S(qs)/qs.length,cyt=S(cs)/cs.length;
   if(Math.min(...qs)<.005||Math.max(...qs)>.995||cyt<.005||cyt>.995)return {status:"boundary_or_transient",q,cyt};
   const prev=seq.slice(0,150),post=seq.slice(150),stable=Math.max(...post.map(v=>v[3]))-Math.min(...post.map(v=>v[3]))<1e-6&&Math.abs(S(prev.map(v=>v[3]))/150-S(post.map(v=>v[3]))/150)<1e-6;
   if(!stable)return {status:"variable_or_transient",q,cyt};
   function lyap(g){
     let z=x.slice(),sum=0,n=0;
     for(let i=0;i<300;i++){let y=rare(x,z,p,g),a=S(y);if(!(a>0))return -Infinity;z=y.map(v=>v/a);if(i>=100){sum+=Math.log(a);n++}}
     return sum/n;
   }
   let lo=-.999,hi=1;while(lyap(hi)<0&&hi<128)hi*=2;
   if(lyap(hi)<0)return {status:"no_bound",q,cyt};
   for(let i=0;i<42;i++){let mid=(lo+hi)/2;if(lyap(mid)>0)hi=mid;else lo=mid}
   let gs=(hi+lo)/2,H=1-q,Ah=p.A/(S(x.map((v,i)=>female(i)?0:v*seeds(i,p)))/H);
   // Nuclear benchmark: same observed q via different adult survival parameters is not equivalent.
   // Sex-only autosomal modifier neutral benchmark, forcing observed q and same mean H seed output.
   // At q not a nuclear equilibrium, do not interpret v68 expression as a valid nuclear invasion threshold.
   if(!(lyap(gs-.0001)<0&&lyap(gs+.0001)>0))throw Error("invasion sign");
   return {status:"stable_interior",q,cyt,realized_A:Ah,gstar:gs,pollinated_self_fraction:p.t,exponent_zero:lyap(0)};
 }
 const all=[];
 for(let A of [1.05,1.15,1.26,1.5,1.75,2,3.77])for(let w of [.3,.5,.7,.9])for(let k of [.1,.3,.5,.7])for(let mode of ["constitutive","silent"])for(let t of [0,.3]){
   let p={A,w,k,mode,t};let a=trajectory(p,[.3,.1,.1,.3,.1,.1]);all.push({...p,...a});
 }
 let states={};for(let r of all)states[r.status]=(states[r.status]||0)+1;
 let stable=all.filter(r=>r.status==="stable_interior"),match=[];
 for(let i=0;i<stable.length;i++)for(let j=i+1;j<stable.length;j++){
   let a=stable[i],b=stable[j];if(a.t!==b.t||Math.abs(a.q-b.q)>.015||Math.abs(Math.log(a.realized_A/b.realized_A))>.12)continue;
   match.push({a,b,diff:Math.abs(a.gstar-b.gstar)})
 }
 match.sort((a,b)=>b.diff-a.diff);
 // Two initial conditions for strongest pair; classify as model witnesses, not fit.
 let witnesses=match.slice(0,5).map(x=>({...x,second_a:trajectory(x.a,[.08,.05,.02,.6,.15,.1]),second_b:trajectory(x.b,[.08,.05,.02,.6,.15,.1])}));
 return {checks,grid_n:all.length,states,stable_n:stable.length,matching_pairs:match.length,witnesses,all};
}
console.log(JSON.stringify(run(),null,2));
