// Deterministic hypothetical population model; no new measured survival or selection.
function run(){
 const eq=(A,sF,sH)=>{const f=1-sF,h=1-sH;return A*h>2*f?(A*h-2*f)/(A*(f+h)-2*f):0};
 const boundary=(A,sF,sH,sM,t)=>{let q=eq(A,sF,sH),h=1-q,T=q*A+h,D=q*(1-sF)+h*(1-sH),O=T-h*t;return 2*h*T*(1-sM)/(D*O)-1};
 const cost=(g,beta,eta)=>{if(g<=0)return 0;let d=(1+g)**(1/beta)-1;return d/(eta+d)};
 const step=(x,A,sF,sH,sM,t,g)=>{
    const surv=[sF,sF,sF,sH,sM,sM];
    const y=x.map((p,i)=>p*surv[i]),D=1-y.reduce((a,b)=>a+b,0);
    const W=x[3]+(1+g)*(x[4]+x[5]),p=(1+g)*(.5*x[4]+x[5])/W;
    const z=Array(6).fill(0);
    for(let i=0;i<6;i++){
        const F=i<3,m=(i%3)/2,seed=x[i]*(F?A:i===3?1-t:0);
        const probs=[(1-m)*(1-p),m*(1-p)+(1-m)*p,m*p];
        for(let j=0;j<3;j++){z[j]+=seed*probs[j]*(F?.5:0);z[j+3]+=seed*probs[j]*(F?.5:1)}
    }
    z[3]+=x[3]*t;
    const total=z.reduce((a,b)=>a+b,0);
    return y.map((v,i)=>v+D*z[i]/total);
 };
 const allele=x=>(x[1]+2*x[2]+x[4]+2*x[5])/2;
 let checks=0,neutralChecks=0;
 for(let A of [1.26,2,3.77,8.0188679245,9.1196741855])for(let sF of [0,.5,.8])for(let sH of [0,.5,.8])for(let sM of [0,.5,.8])for(let t of [0,.25,.5]){
    const q=eq(A,sF,sH),h=1-q,x0=[q,0,0,h,0,0];
    const eqerr=Math.max(...step(x0,A,sF,sH,sM,t,0).map((v,i)=>Math.abs(v-x0[i])));
    if(eqerr>1e-10)throw Error("resident equilibrium failed");
    const gs=boundary(A,sF,sH,sM,t);
    // Positive g+1 multiplier even when threshold is negative.
    const ratio=(g)=>{let eps=1e-10,x=[q,0,0,h-eps,eps,0],start=allele(x);for(let i=0;i<100;i++){x=step(x,A,sF,sH,sM,t,g);if(x.some(v=>!Number.isFinite(v)||v < -1e-12)||Math.abs(x.reduce((a,b)=>a+b,0)-1)>1e-9)throw Error("bad frequency");}return allele(x)/start};
    const lo=ratio(.8*(1+gs)-1),hi=ratio(1.2*(1+gs)-1);
    // Finite-horizon allele counts can initially rise/fall from transient genotype structure.
    // Use spectral invasion condition as primary validation and long-run growth for recurrence.
    const T=q*A+h,D=q*(1-sF)+h*(1-sH),O=T-h*t;
    const lam=(g)=>{if(q===0)return sM+D*(1+g)*O/(2*h*T);const a=sF+D*A/(4*T),b=D*(1+g)*q*A/(4*h*T),c=D*A/(4*T),d=sM+D*(1+g)*(O-q*A/2)/(2*h*T);return (a+d+Math.sqrt((a-d)**2+4*b*c))/2;};
    if(Math.abs(lam(gs)-1)>1e-10||!(lam(.8*(1+gs)-1)<1&&lam(1.2*(1+gs)-1)>1))throw Error("spectral test failed");
    if(!(lo<1&&hi>1))throw Error("recurrence sign failed "+[A,sF,sH,sM,t,q,gs,lo,hi]);
    checks++;
    if(sF===sH&&sH===sM){let q0=A>2?(A-2)/(2*(A-1)):0,T0=q0*A+1-q0,g0=2*(1-q0)*T0/(T0-(1-q0)*t)-1;if(Math.abs(q-q0)>1e-10||Math.abs(gs-g0)>1e-10)throw Error("v68 limit");neutralChecks++}
 }
 let anchor=[["C. palustre",1.26,.193],["C. heterophyllum",3.77,.533]].map(([species,A,qObserved])=>{let relativeLife=(2+qObserved*(A-2))/(A*(1-qObserved));return{species,A,observed_frequency:qObserved,equal_survival_equilibrium:eq(A,0,0),required_F_to_H_expected_adult_census_duration_ratio:relativeLife,scope:"conditional diagnostic: viable achenes assumed adult recruitment, frequency assumed equilibrium, hypothetical dominant nuclear male sterility"}});
 let scenarios=[];
 for(let s of [0,.5,.8])for(let gain of [0,.05,.1]){const A=3.77,t=.25,sF=s,sH=s,sM=s+(1-s)*gain;let gs=boundary(A,sF,sH,sM,t);scenarios.push({A,t,sF,sH,sM,relative_death_reduction:gain,q:eq(A,sF,sH),required_paternity_gain:gs,required_cost_beta1:cost(gs,1,1),required_cost_beta025:cost(gs,.25,1)})}
 return {checks,neutralChecks,anchors:anchor,scenarios,model:"adult vacancy-recruitment recurrence with differential survival; dominant nuclear sex determination; no juvenile lag, stochasticity, CMS or density-dependent pollen limitation"};
}
console.log(JSON.stringify(run(),null,2));
