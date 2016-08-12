load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
begin 
;************************************************ 
; variable and file handling 
;************************************************ 
in = addfile("uwnd.2005.nc","r") ; open netcdf file 
f = addfile("vwnd.2005.nc","r") 
u = short2flt(in->uwnd(152,{200},:,:)) ; pull u off file 
v = short2flt(f->vwnd(152,{200},:,:)) ; pull v off file 
;************************************************ 
; calculate vorticity on a Gaussian Grid 
; scale for visual convenience 
;************************************************ 
scale = 1.e05 
vrt = u ; retain coordinates 
vrt = uv2vrG_Wrap(u,v) * scale 
vrt@long_name = "vorticity" 
vrt@units = "scaled" 
f=2*7.292*1.e-05*sin(phi) 
theta= 300 
pv=g(vrt+f)*d(theta)/dp 
;************************************************ 
; create plot 
;************************************************ 
wks = gsn_open_wks("x11","vort") ; specifies a ps plot 
gsn_define_colormap(wks,"BlWhRe") ; choose colormap 
res = True ; plot mods desired 
res@cnFillOn = True ; turn on color 
res@gsnSpreadColors = True ; use full range of color map 
res@lbLabelAutoStride = True ; nice label bar labels 
res@cnLinesOn = False ; no contour lines 
res@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels 
res@cnMinLevelValF = -5 ; set min contour level 
res@cnMaxLevelValF = 5 ; set max contour level 
res@cnLevelSpacingF = 0.5 ; set contour spacing 
res@lbOrientation = "vertical" 
res@mpMinLatF = 0 ; range to zoom in on 
res@mpMaxLatF = 40. 
res@mpMinLonF = 30. 
res@mpMaxLonF = 120. 
plot = gsn_csm_contour_map_ce(wks,vrt(:,:),res) ; create plot 
end 
