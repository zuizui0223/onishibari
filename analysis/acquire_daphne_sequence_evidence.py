#!/usr/bin/env python3
"""Acquire public, accession-verified Daphne sequences; not an experiment or species tree.
Requires Biopython. Raw source coordinates are NOT exported. Retrieval failures
are explicit and never replaced by invented sequences. Records are dated/hash-bound.
"""
import argparse,csv,hashlib,io,json,re,time,urllib.parse,urllib.request
from datetime import datetime,timezone
from pathlib import Path
from Bio import SeqIO
CP={
 'PV752465':('Daphne pseudomezereum','https://doi.org/10.11110/kjpt.2025.55.3.157'),
 'PV752464':('Daphne kamtschatica','https://doi.org/10.11110/kjpt.2025.55.3.157'),
 'ON244034':('Daphne pseudomezereum var. koreana','https://doi.org/10.1080/23802359.2023.2179356'),
 'MN201546':('Daphne laureola','https://pmc.ncbi.nlm.nih.gov/articles/PMC9209861/'),
 'KY991380':('Daphne kiusiana','https://pmc.ncbi.nlm.nih.gov/articles/PMC9209861/'),
 'MT627481':('Daphne kiusiana var. atrocaulis','https://pmc.ncbi.nlm.nih.gov/articles/PMC9209861/'),
 'BK070978':('Daphne mezereum','https://doi.org/10.1002/2211-5463.70143')}
GENES=('matK','rbcL','rpoB','rpoC1')
BASE='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'

def summarize_record(record,accession,expected):
    if record.id.split('.')[0]!=accession:raise ValueError('accession mismatch')
    organism=record.annotations.get('organism','')
    if expected not in organism:raise ValueError('taxon mismatch: '+organism)
    if expected=='Daphne pseudomezereum' and 'koreana' in organism.lower():
        raise ValueError('koreana cannot substitute for target species')
    seq=str(record.seq).upper()
    if not seq or not set(seq)<=set('ACGTRYSWKMBDHVN'):raise ValueError('invalid DNA')
    genes={}; warnings=[]
    for name in GENES:
        copies=[str(f.extract(record.seq)).upper() for f in record.features
                if f.type=='CDS' and name in f.qualifiers.get('gene',[])]
        if copies and len(set(copies))==1:genes[name]=copies[0]
        else:warnings.append(name+(':missing' if not copies else ':discordant_copies'))
    return dict(accession_version=record.id,organism=organism,length=len(seq),
                sequence_sha256=hashlib.sha256(seq.encode()).hexdigest(),
                genes=genes,warnings=warnings)

def get(url):
    last=None
    for attempt in range(2):
        try:
            req=urllib.request.Request(url,headers={'User-Agent':'onishibari-public-sequence-audit/1.0'})
            with urllib.request.urlopen(req,timeout=35) as response:
                data=response.read(12000000)
                return data,response.geturl(),response.headers.get('Content-Type','')
        except Exception as e:last=e;time.sleep(attempt+1)
    raise RuntimeError(type(last).__name__+': '+str(last)[:250])

def api(endpoint,params):
    return BASE+endpoint+'?'+urllib.parse.urlencode(params)

def fasta(path,label,seq):
    with path.open('a',encoding='utf-8') as fh:
        fh.write('>'+label+'\n'+'\n'.join(seq[i:i+80] for i in range(0,len(seq),80))+'\n')

def main():
    p=argparse.ArgumentParser();p.add_argument('--out',default='outputs/sequence_audit_20261114')
    args=p.parse_args();out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    if (out/'summary.json').exists():raise SystemExit('Use an empty output directory')
    result={'retrieved_utc':datetime.now(timezone.utc).isoformat(),'plastomes':[],
            'its_queries':[],'its_records':[],'resource_attempts':[],
            'limitations':['Not a rooted species tree or a test of mating-system evolution.',
                          'A plastome contains linked genes, not independent nuclear loci.',
                          'Records can represent one individual or cultivation; traits cannot be copied across accessions.']}
    for acc,(expected,doi) in CP.items():
        url=api('efetch.fcgi',dict(db='nuccore',id=acc,rettype='gbwithparts',retmode='text'))
        try:
            raw,final,ctype=get(url); records=list(SeqIO.parse(io.StringIO(raw.decode()),'genbank'))
            if len(records)!=1:raise ValueError('expected single GenBank record')
            rec=records[0];meta=summarize_record(rec,acc,expected)
            meta.update(expected_taxon=expected,source_url=final,source_doi=doi,
                        source_bytes_sha256=hashlib.sha256(raw).hexdigest(),status='verified')
            label=acc+'|'+expected.replace(' ','_')
            fasta(out/'plastomes.fasta',label,str(rec.seq))
            for name,seq in meta['genes'].items():fasta(out/(name+'.fasta'),label,seq)
            meta['gene_lengths']={k:len(v) for k,v in meta.pop('genes').items()}
            result['plastomes'].append(meta)
            print('PLASTOME='+json.dumps(meta,ensure_ascii=False),flush=True)
        except Exception as e:
            error=dict(accession=acc,expected_taxon=expected,source_url=url,status='failed',error=str(e))
            result['plastomes'].append(error);print('PLASTOME='+json.dumps(error),flush=True)
        time.sleep(.45)
    seen=set()
    for taxon in ('Daphne pseudomezereum','Daphne kamtschatica','Daphne jezoensis','Daphne arisanensis'):
        term=f'"{taxon}"[Organism] AND ("internal transcribed spacer"[Title] OR ITS1[Title] OR ITS2[Title]) AND 200:5000[SLEN]'
        url=api('esearch.fcgi',dict(db='nuccore',term=term,retmode='json',retmax=40))
        q={'query_taxon':taxon,'query':term,'query_url':url}
        try:
            raw,_,_=get(url);ans=json.loads(raw)['esearchresult'];ids=ans['idlist']
            q.update(count=int(ans['count']),returned_ids=len(ids),truncated=int(ans['count'])>len(ids),status='queried')
            if ids:
                time.sleep(.45)
                furl=api('efetch.fcgi',dict(db='nuccore',id=','.join(ids),rettype='gb',retmode='text'))
                data,final,_=get(furl)
                for rec in SeqIO.parse(io.StringIO(data.decode()),'genbank'):
                    if rec.id in seen:continue
                    seen.add(rec.id);organism=rec.annotations.get('organism','')
                    seq=str(rec.seq).upper()
                    if not organism.startswith('Daphne ') or not set(seq)<=set('ACGTRYSWKMBDHVN'):continue
                    meta=dict(accession_version=rec.id,organism=organism,length=len(seq),
                              description=rec.description,source_query_taxon=taxon,
                              source_url=final,sequence_sha256=hashlib.sha256(seq.encode()).hexdigest())
                    result['its_records'].append(meta)
                    fasta(out/'its_candidates.fasta',rec.id+'|'+organism.replace(' ','_'),seq)
            result['its_queries'].append(q);print('ITS_QUERY='+json.dumps(q),flush=True)
        except Exception as e:q.update(status='failed',error=str(e));result['its_queries'].append(q);print('ITS_QUERY='+json.dumps(q),flush=True)
        time.sleep(.45)
    resources={
      'treebase_S10851':'https://purl.org/phylo/treebase/phylows/study/TB2:S10851?format=nexus',
      'han2025_pdf':'https://e-kjpt.org/upload/pdf/kjpt-2025-55-3-157.pdf',
      'alonso2011_pdf':'https://www.plant-animal.es/pdfs/Alonso%26Herrera_Evolution_2011.pdf',
      'ramos2024_full':'https://onlinelibrary.wiley.com/doi/full/10.1111/jse.13018'}
    for label,url in resources.items():
        m=dict(label=label,url=url)
        try:
            data,final,ctype=get(url);m.update(final_url=final,bytes=len(data),sha256=hashlib.sha256(data).hexdigest(),content_type=ctype)
            is_tree=data.lstrip().upper().startswith(b'#NEXUS')
            is_pdf=data.lstrip().startswith(b'%PDF')
            if is_tree or is_pdf:
                name=label+('.nex' if is_tree else '.pdf');(out/name).write_bytes(data)
                m.update(status='downloaded',file=name)
            else:
                m.update(status='not_phylogenetic_file',preview=data[:100].decode(errors='replace'))
        except Exception as e:m.update(status='failed',error=str(e))
        result['resource_attempts'].append(m);print('RESOURCE='+json.dumps(m),flush=True)
    result['verified_plastomes']=sum(x['status']=='verified' for x in result['plastomes'])
    result['its_record_count']=len(result['its_records'])
    (out/'summary.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    with (out/'its_inventory.csv').open('w',newline='') as fh:
        fields=['accession_version','organism','length','description','source_query_taxon','source_url','sequence_sha256']
        w=csv.DictWriter(fh,fieldnames=fields);w.writeheader();w.writerows(result['its_records'])
    print('AUDIT_SUMMARY='+json.dumps({k:result[k] for k in ('verified_plastomes','its_record_count')}),flush=True)
    if result['verified_plastomes']==0:raise SystemExit('No verified plastomes; no successful sequence acquisition claim allowed')
if __name__=='__main__':main()
