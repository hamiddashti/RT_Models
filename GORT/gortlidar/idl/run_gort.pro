pro run_gort

 tek_color
  ; Setting plotting
  set_plot,'ps'
  device, file= 'comp.eps', /color,/TIMES,/BOLD,/portrait, ysize=28, YOFFSET=0, xsize=21, xOFFSET=0.0
  !p.multi=[0,2,3]

stnamein=['noG.BAC2_SW','noG.BAC2_SE'] ; noGaussian smooth

stnameout=['G.BAC2_SW','G.BAC2_SE'] ; Gaussian smoothed 

nsite=2
indir='../output_noG/'
outdir='../output_G/'

for i=0, nsite-1 do begin

  ;run guassian convolution
  temp_str=''
  gn_level=0
  gout = {height:0.0,fp:0.0,efp:0.0,pgap_noc:0.0,pgap_c:0.0,dpdz_noc:0.0,dpdz_c:0.0,wvfm_noc:0.0,wvfm_c:0.0}

;read gort output file 
  close,2
  openr, 2, indir+stnamein(i)
  readf, 2, gn_level,dz
  n_ext=10./dz  
  n_tot=gn_level+n_ext

  ght=fltarr(n_tot)
  wvfm_noc=fltarr(n_tot)
  wvfm_c=wvfm_noc
  gwvfm_noc=wvfm_noc
  gwvfm_c=wvfm_noc

  ght(0:n_ext-1)=findgen(n_ext)*dz-10 ;extend height to below ground 10m deep 
  wvfm_noc(0:n_ext-1)=0             ;extend wvfm to below ground 10m deep
  wvfm_c(0:n_ext-1)=0             ;extend wvfm to below ground 10m deep  

  readf, 2, temp_str
  readf, 2, temp_str
  readf, 2, temp_str
  for ilevel=0, gn_level-1 do begin
;    print,ilevel
      readf, 2, gout
      ght[ilevel+n_ext]=gout.height
      wvfm_noc[ilevel+n_ext]=gout.wvfm_noc
      wvfm_c[ilevel+n_ext]=gout.wvfm_c
  endfor
  close, 2

  width = fix(0.6/dz)+1 ;LVIS gaussian width=10ns(10ns*0.15m=1.5m)/2.35482=0.63m, witdth #=0.63/dz
  gwvfm_noc=gsmooth(wvfm_noc,width)
  gwvfm_c=gsmooth(wvfm_c,width)
 
  if (i eq 0) then begin
      plot, gwvfm_noc, ght,charsize=2
      oplot, gwvfm_c, ght,color=2
      print, 'isite=',i
      endif else begin  
          plot, gwvfm_noc, ght, charsize=2
          oplot, gwvfm_c, ght, color=2
          print, 'isite=',i 
          endelse
 
  save, filename=outdir+stnameout(i)+'.sav',ght,gwvfm_noc,gwvfm_c
endfor

DEVICE,/CLOSE 
stop

end
