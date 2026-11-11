void cmpv2(){ const char* v[3]={"kfp","raw","kappa"}; const char* h[6]={"hKs_qa","hKs_loose","hLa_qa","hLa_loose","hLaBar_qa","hLaBar_loose"};
 int nf=0, bad=0, miss=0; TSystemDirectory d("d","v0mass_kappa_v2"); TList* L=d.GetListOfFiles();
 for (auto o: *L) { TString n=o->GetName(); if(!n.EndsWith("_kappa.root")) continue; TString b=n; b.ReplaceAll("_kappa.root","");
  for(int i=0;i<3;++i){ TFile x(Form("v0mass_kappa_v2/%s_%s.root",b.Data(),v[i])); TFile y(Form("v0mass_kappa/%s_%s.root",b.Data(),v[i])); if(y.IsZombie()){++miss;continue;}
   for(int j=0;j<6;++j){ TH2* p=(TH2*)x.Get(h[j]); TH2* q=(TH2*)y.Get(h[j]); for(int k=0;k<(p->GetNbinsX()+2)*(p->GetNbinsY()+2);++k) if(p->GetBinContent(k)!=q->GetBinContent(k)){++bad; printf("DIFF %s %s %s\n",b.Data(),v[i],h[j]); break;} } }
  ++nf; }
 printf("REPROV2 segments %d differing-histograms %d missing-old %d\n", nf, bad, miss); }
