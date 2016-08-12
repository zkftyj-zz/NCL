load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
	path		= "../download/"

;=====================
	iyear = 2004

	fils 		= systemfunc ("csh -c 'cd "+path+" ; ls *"+iyear+"*.he5'")
	nf		= dimsizes(fils)
	f    	= addfile(path+fils(0),"r")
	print(nf)

	ft		= addfile("/home/u1/zk794/lonestar/research/data/Aura_MLS/ncar/data/mls_h2o_daily_v3.nc","r")
	lat		= ft->lat
	lon		= ft->lon
	nlat		= dimsizes(lat)
	nlon		= dimsizes(lon)
	
	indlev		= ind(f->Pressure_H2O.ge.0.002 .and. f->Pressure_H2O.lt.316.)
	lev		= f->Pressure_H2O(indlev)
	nlev		= dimsizes(lev)

	do i = 0,nf-1
		if (isfilepresent(fils(i)+".nc"))
			continue
		else
		xx 			= new((/nlev,nlat,nlon/),float,-999)
		f    = addfile(path+fils(i),"r")

		do	ilev	= 0,nlev-1
			x	= new((/nlat,nlon/),float,-999)
			N	= new((/nlat,nlon/),float,-999)
			x 	= 0.
			N 	= 0.

			vlev	= lev(ilev)
			ilevind	= ind(f->Pressure_H2O.eq.vlev)

			lat0	= f->Latitude_H2O
			lon0	= where(f->Longitude_H2O.lt.0,f->Longitude_H2O+360.,f->Longitude_H2O)
			x0	= f->L2gpValue_H2O(:,ilevind)
			
			cv	= f->Convergence_H2O
			qu	= f->Quality_H2O
			pr	= f->L2gpPrecision_H2O(:,ilevind)
			st	= f->Status_H2O
			
			indx	= ind((cv.lt.2.0) .and. (qu.gt.1.45) .and. (mod(st,2).eq.0) .and. (pr.gt.0).and. (.not.ismissing(x0)).and. (.not.ismissing(lat0)).and. (.not.ismissing(lon0)))
			nindx	= dimsizes(indx)

			if (nindx.gt.2)
				do	iindx = 0,nindx-1
				latind	= lat0(indx(iindx))
				lonind	= lon0(indx(iindx))

				xlatind	= closest_val(latind,lat)
				xlonind	= closest_val(lonind,lon)
	
				if	(xlonind.gt.nlon-1)
					xlonind        = xlonind-1
				end if

				x(xlatind,xlonind)	= x(xlatind,xlonind) + x0(indx(iindx))
				N(xlatind,xlonind) = N(xlatind,xlonind) + 1.
				delete([/latind,lonind,xlatind,xlonind/])
				end do
			end if

			None		= ndtooned(N)
			xone		= ndtooned(x)
			indN		= ind(None.gt.0)
			indN0		= ind(None.eq.0)
			xone(indN0)	= -999
			if (dimsizes(indN).gt.1)
			xone(indN)	= xone(indN)/None(indN)
			end if
			N		= onedtond(None,(/nlat,nlon/))
			xx(ilev,:,:)		= onedtond(xone,(/nlat,nlon/))

			delete([/x,N,lat0,lon0,x0,ilevind,cv,qu,pr,st,indx,nindx,None,xone,indN,indN0/])
		end do

		xx!0 		= "lev"
		xx!1 		= "lat"
		xx!2 		= "lon"
		xx&lev	= f->Pressure_H2O(indlev)
		xx&lat	= lat
		xx&lon	= lon
		printVarSummary(xx)

		fon		= fils(i)+".nc"
		system("rm -rf "+fon)
		fo		= addfile(fon,"c")
	
		fo->H2O		= xx
		xtime 		= dim_avg(f->Time_H2O)
		xtime@units 	= "seconds since 1993-1-1 00:00:00"
		fo->time 	= xtime
		print(fon+"-----complete!-------")
		delete([/xx,xtime/])
		end if
	end do

end 
