
f2data = READ_CSV('csv\bO3.csv',N_TABLE_HEADER=1)


; Variables to be taken from model.inp
ICE_AREA = 1.0e-20 ; cm2
ICE_THICK = 1.0e-4 ; cm
ICE_DENSITY = 5.7e22 ;molecules/cm3 
flux = 2.33E14

; Calculate initial O2 
initialO2 = ICE_DENSITY*ICE_AREA*ICE_THICK ; # of initial O2 molecules
initial_oxygen = ICE_DENSITY*ICE_AREA*ICE_THICK ; # of initial O atoms (in O2 molecules)

set_plot, 'WIN'

expX = [4.709604,$
13.496563,$
50.07993,$
146.55681,$
465.09692,$
2581.1648]

expY = [0.055385794,$
0.15238622,$
0.52023816,$
1.062651,$
1.4836915,$
1.6182984]

expY = expY*15
expX = expX*flux

f2data.field1 = f2data.field1
;f2data.field2 = ((f2data.field2*ICE_THICK*ICE_AREA)/initial_oxygen)*100
f2data.field2 = (f2data.field2/ice_density)*100


;DEVICE, FILENAME='f2.eps', DECOMPOSED=0, /ENCAPSULATED

plot, f2data.field1,$
      f2data.field2,$
      /XLOG,$
      /YLOG,$
     XRANGE=[1e15,1e18],$
      YRANGE=[0.8,50],$
      xtitle='Fluence (photons/cm!E2!N)',$
      ytitle='[O!D3!N]/[initial O!D2!N] x 100 %',$
;      title='O!D3!N Production During O!D2!N Irradiation',$
      ystyle=1,$
      charsize=1.8,$
      charthick=3,$
      xthick=3,$
      ythick=3,$
      /NODATA

loadct, 40

      
oplot, expX,$
       expY,$
       PSYM=-5,$
       COLOR=250,$
       thick=5,$
       symsize=1.5
       
oplot, f2data.field1,$
       f2data.field2,$
       COLOR=250,$
       thick=5

legend,['Experiment','Model'],psym=[-5,0],number=1,/bottom,/right,colors=[250,40],box=0,charsize=1.8,charthick=3,thick=5,pspacing=1          ; plot two symbols, not one
      
;device, /close

END
