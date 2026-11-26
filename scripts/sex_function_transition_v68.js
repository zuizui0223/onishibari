// Node.js; no external packages. Counterfactual only. All inputs below are published point summaries; no individual data or new observed fitness.
const review=[
  [
    "Beta vulgaris",
    0.93,
    "NS",
    "NC"
  ],
  [
    "Gingidia flabellata",
    1.49,
    "NFL*FS",
    "unknown"
  ],
  [
    "Echium vulgare",
    1.42,
    "NS",
    "NC"
  ],
  [
    "Phacelia linearis",
    1.68,
    "NS",
    "N"
  ],
  [
    "Opuntia quimilo",
    0.88,
    "NS",
    "unknown"
  ],
  [
    "Lobelia siphilitica",
    1.31,
    "NS",
    "NC"
  ],
  [
    "Dianthus sylvestris",
    1.62,
    "NFL*NSF",
    "unknown"
  ],
  [
    "Silene italica",
    1.22,
    "NS",
    "NC"
  ],
  [
    "Silene vulgaris",
    2.17,
    "FS*NSF",
    "NC"
  ],
  [
    "Wurmbea biglandulosa",
    1.77,
    "NS",
    "unknown"
  ],
  [
    "Cucurbita foetidissima",
    1.5,
    "NS",
    "N"
  ],
  [
    "Geranium maculatum",
    1.6,
    "NS",
    "unknown"
  ],
  [
    "Geranium sylvaticum",
    1.68,
    "NS",
    "unknown"
  ],
  [
    "Nemophila menziesii",
    1.67,
    "NS",
    "NC"
  ],
  [
    "Lycopus maackianus",
    37.73,
    "NS",
    "unknown"
  ],
  [
    "Thymus vulgaris",
    2.5,
    "NS",
    "NC"
  ],
  [
    "Chionographis japonica",
    1.4,
    "NS",
    "unknown"
  ],
  [
    "Hebe subalpina",
    3.44,
    "NFL*FS",
    "unknown"
  ],
  [
    "Plantago coronopus",
    1.72,
    "NS#",
    "NC"
  ],
  [
    "Plantago lanceolata",
    1.38,
    "NS",
    "NC"
  ],
  [
    "Plantago maritima",
    1.01,
    "NS",
    "NC"
  ],
  [
    "Chionochloa bromoides",
    0.91,
    "NFL*NSF",
    "unknown"
  ],
  [
    "Fragaria virginiana",
    6.68,
    "NFL*FS",
    "N"
  ],
  [
    "Prunus mahaleb",
    1.4,
    "NF",
    "unknown"
  ],
  [
    "Saxifraga granulata",
    0.72,
    "NFL*NSF",
    "NC"
  ],
  [
    "Daphne laureola",
    0.65,
    "NS",
    "unknown"
  ],
  [
    "Gnidia wikstroemiana",
    3.58,
    "NFL*FS",
    "unknown"
  ],
  [
    "Kallstroemia grandiflora",
    1.29,
    "NS",
    "unknown"
  ]
];
const anchors=[
  {
    "species": "Cirsium palustre",
    "A": 1.26,
    "female_frequency": 0.193,
    "endpoint": "viable achenes per plant",
    "frequency_scope": "mean over 285 populations",
    "doi": "10.1002/ece3.74093"
  },
  {
    "species": "Cirsium heterophyllum",
    "A": 3.77,
    "female_frequency": 0.533,
    "endpoint": "viable achenes per plant",
    "frequency_scope": "mean over 126 populations",
    "doi": "10.1002/ece3.74093"
  },
  {
    "species": "Daphne jezoensis 2018",
    "A": 8.018867924528301,
    "female_frequency": 0.375,
    "endpoint": "natural fruit per flowering genet",
    "frequency_scope": "strict F/H flowering genets with valid fruit records, NOT population sex ratio",
    "doi": "10.1002/ajb2.16239"
  },
  {
    "species": "Daphne jezoensis 2020",
    "A": 9.11967418546366,
    "female_frequency": 0.4453125,
    "endpoint": "natural fruit per flowering genet",
    "frequency_scope": "strict F/H flowering genets with valid fruit records, NOT population sex ratio",
    "doi": "10.1002/ajb2.16239"
  }
];
function run(review,anchors){
  const q=A=>A>2?(A-2)/(2*(A-1)):0;
  const gstar=(A,t)=>{const f=q(A),h=1-f,T=f*A+h,O=T-h*t;return 2*h*T/O-1;};
  const cstar=(g,beta,eta)=>{const d=(1+g)**(1/beta)-1;return d/(eta+d);};
  // State [F BB,F Bb,F bb,H BB,H Bb,H bb], F is Aa; H aa.
  // b is an unlinked dominant pistil-loss modifier in aa plants only.
  // t is the fraction of viable H BB maternal recruits arising by selfing.
  function step(x,A,t,g){
    const h=x[3]+(1+g)*(x[4]+x[5]);
    const p=(1+g)*(.5*x[4]+x[5])/h;
    const y=Array(6).fill(0);
    for(let i=0;i<6;i++){
      const female=i<3, m=(i%3)/2;
      const seed=x[i]*(female?A:i===3?1-t:0);
      const z=[(1-m)*(1-p),m*(1-p)+(1-m)*p,m*p];
      for(let j=0;j<3;j++){
        y[j]+=seed*z[j]*(female?.5:0);
        y[j+3]+=seed*z[j]*(female?.5:1);
      }
    }
    y[3]+=x[3]*t;
    const total=y.reduce((a,b)=>a+b,0);
    return y.map(v=>v/total);
  }
  const allele=x=>(x[1]+2*x[2]+x[4]+2*x[5])/2;
  function check(A,t,g){
    let f=q(A),h=1-f, eps=1e-9;
    let x=[f,0,0,h-eps,eps,0],start=allele(x);
    for(let i=0;i<100;i++)x=step(x,A,t,g);
    if(Math.abs(x.reduce((a,b)=>a+b,0)-1)>1e-10||x.some(v=>v<0))throw Error("invalid population");
    return allele(x)/start;
  }
  const rows=[...review.map(([species,A,endpoint,genetics])=>({species,A,endpoint,genetics,source:"Dufay & Billard table 6; DOI 10.1093/aob/mcr062"})),...anchors];
  const results=rows.map(r=>({...r,q_nuclear_counterfactual:q(r.A),gstar_t0:gstar(r.A,0),gstar_t50:gstar(r.A,.5),
    cost_t0_beta1:cstar(gstar(r.A,0),1,1),cost_t0_beta025:cstar(gstar(r.A,0),.25,1),
    cost_t50_beta1:cstar(gstar(r.A,.5),1,1),cost_t50_beta025:cstar(gstar(r.A,.5),.25,1)}));
  let checks=0;
  for(const r of rows)for(const t of [0,.25,.5,.75]){
    const f=q(r.A),x=[f,0,0,1-f,0,0],y=step(x,r.A,t,0);
    if(Math.abs(y[0]-f)>1e-12)throw Error("resident equilibrium failed");
    const g=gstar(r.A,t),below=check(r.A,t,g*.5),above=check(r.A,t,g*1.5);
    if(!(below<1&&above>1))throw Error("invasion sign failed "+r.species+" "+t+" "+below+" "+above);
    checks++;
  }
  return {review_n:review.length,review_A_gt2:review.filter(r=>r[1]>2).length,
    source_reported_nuclear:results.filter(r=>r.genetics==="N").map(r=>({species:r.species,A:r.A,q:r.q_nuclear_counterfactual})),
    checks,model:"counterfactual dominant nuclear male sterility + unlinked dominant pistil-loss modifier; fixed relative viable maternal recruitment; random outcrossing; nonoverlapping generations; t is realized viable selfed fraction, not raw selfing rate",
    rows:results};
}
console.log(JSON.stringify(run(review,anchors),null,2));
