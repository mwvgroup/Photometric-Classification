templatetype 91bgnugent
maxs 2.6
daymax 85
forcefit and fredshift=sne.redshift
wanted_filters=['u','g','r','i','z']+'_2.5m'


CALL
newsub,snname=sne[w].name,templatetype='91bgnlmod',/nosnake,daymax=85.,verbose=1,/forcef,$
fittype='simple',wanted_filters=['u','g','r','i','z']+'_2.5m',namedir='sdss-91bg',/ffit,$ 
fredshift=sne[w].redshift,survey='sdss',maxs=2.6                
