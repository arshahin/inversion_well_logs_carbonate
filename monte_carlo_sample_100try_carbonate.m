function sample=monte_carlo_sample_100try_carbonate(phimin,dphi,phit,phiiw,phivw)
nphi=round(abs(phit-phimin)./dphi);
ntry=1;
ip=2.*(rand-0.5).*nphi;
sample=phimin+ip.*dphi;
while (sample < phimin || sample > phit || sample < phiiw || sample > (phit-phivw)) 
    ntry=ntry+1;
    ip=2.*(rand-0.5).*nphi;
    sample=phimin+ip.*dphi;
    if(ntry>1000) 
    break
    end
end
end
    







