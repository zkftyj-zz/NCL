load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

; This function is written to derive the partial correlation between X,Y,and Z(n,:)
; X & Y are one-dimension matrix with "time" ntime, and Z is two-dimension: Z(n,ntime)
; Written by Kai Zhang

function partial(n[1],X[*]:float,Y[*]:float,Z[*][*]:float)
local rxy,rxz,ryz,partial1,Rxy,Rxz,Ryz,partialn
begin
	if(n.eq.1) then
		rxy		= escorc(X,Y)
		rxz		= escorc(X,Z)
		ryz		= escorc(Y,Z)
		sqrtpar	= sqrt((1.-rxz*rxz)*(1.-ryz*ryz))
		if	(.not.ismissing(sqrtpar).and.sqrtpar.eq.0.) then
		partial1	= default_fillvalue("float")
		else if	(ismissing(rxy).or.ismissing(rxz).or.ismissing(ryz).or.ismissing(sqrtpar)) then
		partial1	= default_fillvalue("float")
		else
		partial1	= (rxy-rxz*ryz)/sqrtpar
		end if
		end if
		return(partial1)
	else
		Rxy		= partial(n-1,X,Y,Z(0:n-2,:))
		Rxz		= partial(n-1,X,Z(n-1,:),Z(0:n-2,:))
		Ryz		= partial(n-1,Y,Z(n-1,:),Z(0:n-2,:))
		if	(ismissing(Rxy).or.ismissing(Rxz).or.ismissing(Ryz))	then
		partialn	= default_fillvalue("float")
		else
		partialn	= (Rxy-Rxz*Ryz)/sqrt((1.-Rxz*Rxz)*(1.-Ryz*Ryz))
		end if
		return(partialn)
	end if
end
