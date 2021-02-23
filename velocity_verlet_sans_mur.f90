	Program velocity_verlet
	Implicit none

	integer:: n,pot,i
	integer:: iseed	
	real(8),allocatable::x(:), y(:), z(:) 
	real(8),allocatable::vx(:), vy(:), vz(:) 
	real(8)::v0
!	real(8),allocatable::ax(:), ay(:), az(:) 
	real(8),allocatable::fx(:), fy(:), fz(:) 
	real(8),allocatable::vmod(:)
	real(8)::velcmx,velcmy,velcmz
	real(8)::t,tmax
	real(8)::dt,densidade,raiodecorte	
	real(8),parameter::Epsi=1.0d0,Sigma=1.0d0
	real(8):: ktot,utot
	real(8):: pi,tupi	
	real(8):: modulotempo
	integer,parameter::tentativa_max=2000
	real(8)::L,dmin,dmin2
	real(8)::tterm
	integer::numerociclos
	integer::cicloatual
	integer::numerodequadros
	integer::numeropassosporquadro
!	integer::numerolados
!	character(len = 3)::tipopotencial

!	print*,'Qual o tipo de potencial? '
!	print*,'(1) para Lennard-Jones (2) WCA'
	read*,pot
!	print*,'Qual o intervalo de tempo?'
	read*,dt
!	print*,'Durante quanto tempo?'
	read*,tmax
!	print*,'quanto tempo para relaxar?'
	read*,tterm
!	print*,'Quantas Particulas?'
	read*,N
!	print*,'Entre um inteiro negativo para seed'
	read*,iseed
!	print*,' qual a densidade dessa bagaça?'
	read*,densidade
!	print*, 'Qual a separacao mínima entre partículas?'
	read*,dmin
!	print*, 'Qual o modulo da vel. inicial das partículas?'
	read*,v0
!	print*, 'Queres comprimir o tempo de vossa simulação? se sim por quantas vezes?
	read*,numerodequadros
	
	L=(n/densidade)**(1.0d0/3.0d0)
	if (pot	== 1) then
		raiodecorte = 2.5d0
		print*,'Potencial LJ de um sistema de ',n,'particulas em um volume de lado',L,','
	else
		raiodecorte = 2.0d0**(1.0d0/6.0d0)		
		print*,'Potencial WCA de um sistema de ',n,'particulas em um volume de ',L,','
	end if
		print*,'tempo entre posicoes dt =',dt,', a parede caiu'			
	

	allocate(x(n), y(n), z(n))
	x = 0.0d0
	y = 0.0d0
	z = 0.0d0

	allocate(fx(n), fy(n), fz(n))
	fx = 0.0d0
	fy = 0.0d0
	fz = 0.0d0

	allocate(vx(n), vy(n), vz(n), vmod(n))
	vx = 0.0d0
	vy = 0.0d0
	vz = 0.0d0
	vmod = 0.0d0
	

	pi = 4.0d0 * atan(1.0d0)
	tupi = 2.0d0 * pi

	open(unit=10,file="saida.dat")
	open(unit=11,file="filme.xyz")
    open(unit=13,file='plotvelcm.plt')
		write(13,*)'set terminal postscript eps enhanced color font "Helvetica,10"'
		write(13,*)'set output "centrodemassa.eps"'
		write(13,*)'set xlabel "tempo(s)"'
		write(13,*)'set ylabel "Velocidade(m/s)"'
!		write(13,*)'set zlabel "PSIxy"'
		!write(13,*)'set pm3d'
!		write(13,*)'set view 30,30'
!		write(13,*)'set xrange [0:1]'
!		write(13,*)'set yrange [0:1]'		
!		write(13,*)'set dgrid3d 50 50'
		!write(13,*)'hidden3d'
		write(13,*)'m = "saida.dat"'
		write(13,*)'set encoding utf8'
		write(13,*)'set terminal x11 0'
		!write(13,*)'nokey'
		!write(13,*)'grid'
		!write(13,*)'set xtics ("100" 100, "200" 200, "300" 300)'
		write(13,*)'set title "Gráfico Velocidade centro de massa"'
		!write(13,*)'set xtics ("100" 0, "200" 1, "300" 2)'
		!write(13,*)'set ytics ("100" 0, "200" 1, "300" 2)'
		!write(13,*)'splot "dados.dat" u 1:2:3 with lines lc rgb "red"' !'plot m using 1:2 with linespoints lc rgb "blue"' !, m using 3:1 with linespoints lc rgb "red" '
		write(13,*)'plot "saida.dat" u 1:2 with lines lc rgb "red", \'
		write(13,*)'	"saida.dat" u 1:3 with lines lc rgb "green" ,\' 
!		write(13,*)'	"saida.dat" u 1:6 with lines lc rgb "blue" '
	close(13)

!-------------------------Programa Principal-----------------------
	densidade = n/L**3
	
	call init ()
	call force (1)
	t = 0.0d0 !zera o tempo
	print*, utot
	write(10,'(9(ES13.5E3,1x))')t,utot/real(n),ktot/real(n),velcmx,velcmy,velcmz,vmod(10),vmod(50),vmod(70)
	write(11,'(1x,i6)')n
	write(11,'(A)')'filmeee'
	do i = 1,n
		write(11,'(1x,A1,1x,ES12.5,1x,ES12.5,1x,ES12.5)')'O',x(i),y(i),z(i)
	end do
	cicloatual = 0
	numerociclos=int(tmax/dt)
	numeropassosporquadro=int(numerociclos/numerodequadros)
	do  while (t < tmax)
	!	modulotempo = modulo(cicloatual,numerociclos*compressaotemporal)
		call integrate (1)
		cicloatual = cicloatual + 1
		if (modulo(cicloatual,numeropassosporquadro) == 0) then
			write(11,'(1x,i6)')n
			write(11,'(A)')'filmeee'
			do i = 1,n
				write(11,'(1x,A1,1x,ES12.5,1x,ES12.5,1x,ES12.5)')'O',x(i),y(i),z(i)
			end do
		end if
		call force (1)
		call integrate (2)
		if (t > tterm) then
			call integrate (3)
			call force (2)
			write(10,'(9(ES12.5,1x))')t,utot/real(n),ktot/real(n),velcmx,velcmy,velcmz,vmod(10),vmod(50),vmod(70)
		end if
		t = t + dt
!		call sample ()
	end do
    call system('gnuplot -p plotvelcm.plt' ) !chamar gnuplot usando script plotvelcm.plt
	call system('epstopdf centrodemassa.eps') !converter pdf
	call system('centrodemassa.pdf') !abrir pdf
	close(unit=11)
!------------------------Programa Principal------------------------

	close(unit=10)
	CONTAINS

	Subroutine init ()
	real(8) :: modraio
	real(8) :: teta,phi
	real(8) :: dx, dy ,dz
	integer :: i,j
	logical :: overlap
	integer :: tentativa	
	
	tentativa = 1
	i = 1
	do while (i.le.n)
!		j=i+1	
		x(i) =  ran2(iseed)*L
		!print*,iseed
		y(i) = ran2(iseed)*L
		!print*,iseed
		z(i) = ran2(iseed)*L
		j = 1
		overlap = .false.
		do while (j < i  .and. (.not. overlap))
			dx=x(i)-x(j)
			dy=y(i)-y(j)
			dz=z(i)-z(j)
			dx = dx - L * ANINT(dx / L)
			dy = dy - L * ANINT(dy / L)
			dz = dz - L * ANINT(dz / L)
			modraio=SQRT( dx * dx + dy * dy + dz * dz)
					
			if (modraio < dmin) then		
				overlap = .true.
			end if
			j=j+1
		end do
		if (overlap) then
			tentativa = tentativa +1
			if (tentativa > tentativa_max) then
				print*,' HERESIA!!!Tentativas maximas excedidas!!'
				stop
			endif
		else
			tentativa = 1
			i = i + 1
		endif
	end do
!	velcmx = 0.0d0
!	velcmy = 0.0d0
!	velcmz = 0.0d0
	vx(i) = 0.0d0
	vy(i) = 0.0d0
	vz(i) = 0.0d0
	ktot = 0.0d0
	vmod = 0.0d0
	do i = 1, n
		teta = ran2(iseed)*pi
		phi  = ran2(iseed)*tupi	
		vx(i) = v0*sin(teta)*cos(phi)
		vy(i) = v0*sin(teta)*sin(phi)
		vz(i) = v0*cos(phi)
		vmod(i) = SQRT(vx(i)*vx(i) + vy(i)*vy(i) + vz(i)*vz(i))
		ktot = ktot + 0.5d0 * (vx(i)*vx(i) + vy(i)*vy(i) + vz(i)*vz(i))
		velcmx = velcmx + vx(i)
		velcmy = velcmy + vy(i)
		velcmz = velcmz + vz(i)
	end do 
	velcmx = velcmx/real(n)
	velcmy = velcmy/real(n)
	velcmz = velcmz/real(n)

!	do i=1,n   		!isso é um truque barato, em sistema fechado(parede)
!	vx(i) = vx(i) - velcmx  !
!	vy(i) = vy(i) - velcmy
!	vz(i) = vz(i) - velcmz
	
!	velcmx = velcmx + vx(i)
!	velcmy = velcmy + vy(i)
!	velcmz = velcmz + vz(i)

!	end do
	End Subroutine init

!-------------------------Switch entre integracoes-----------------



	subroutine integrate (switch)
	integer :: switch
	integer :: i

		if (switch == 1) then
			do i = 1, N
!					vx(i) = vx(i) + 0.5 * ax(i) * dt vamos fazer diferente
!					vy(i) = vy(i) + 0.5 * ay(i) * dt usaremos coordenadas polares
!					vz(i) = vz(i) + 0.5 * az(i) * dt como segue abaixo:
					vx(i) = vx(i) + 0.5d0 * fx(i) * dt
					vy(i) = vy(i) + 0.5d0 * fy(i) * dt
					vz(i) = vz(i) + 0.5d0 * fz(i) * dt
					x(i) = x(i) + vx(i) * dt
					y(i) = y(i) + vy(i) * dt
					z(i) = z(i) + vz(i) * dt
					if (x(i) > L) then
						x(i) = x(i) - L
					else if (x(i) < 0.0d0) then
						x(i) = x(i) + L
 					end if
					if (y(i) > L) then
						y(i) = y(i) - L
					else if (y(i) < 0.0d0) then
						y(i) = y(i) + L
 					end if
					if (z(i) > L) then
						z(i) = z(i) - L
					else if (z(i) < 0.0d0) then
						z(i) = z(i) + L
 					end if
			end do
		else if (switch == 2) then
			do i = 1, N
					vx(i) = vx(i) + 0.5d0 * fx(i) * dt
					vy(i) = vy(i) + 0.5d0 * fy(i) * dt
					vz(i) = vz(i) + 0.5d0 * fz(i) * dt	
			end do				
		else if (switch == 3) then
			velcmx = 0.0d0
			velcmy = 0.0d0
			velcmz = 0.0d0
			ktot = 0.0d0	
			do i = 1, N
					vmod(i) = SQRT(vx(i)*vx(i) + vy(i)*vy(i) + vz(i)*vz(i))
					ktot = ktot + 0.5d0 * (vmod(i)*vmod(i))
					velcmx = velcmx + vx(i)
					velcmy = velcmy + vy(i)
					velcmz = velcmz + vz(i)
			end do				
			velcmx = velcmx/real(n)
			velcmy = velcmy/real(n)
			velcmz = velcmz/real(n)
		else
			print*, "deu pau"
			stop
		end if
	end subroutine integrate


	subroutine force (switch)
		integer::switch
		integer::i,j
		real(8)::modraio,factor
		real(8)::dr(3)		
		if (switch == 1) then
		fx = 0.0d0
		fy = 0.0d0
		fz = 0.0d0
		i = 1
		
		do while (i.le.n-1)
			j = i + 1
			do while (j.le.n)	
				dr(1) = x(i) - x(j)
				dr(2) = y(i) - y(j)
				dr(3) = z(i) - z(j)
				dr(1) = dr(1) - L * ANINT(dr(1) / L)
				dr(2) = dr(2) - L * ANINT(dr(2) / L)
				dr(3) = dr(3) - L * ANINT(dr(3) / L)
				modraio=SQRT( dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3))
				if (modraio.le.raiodecorte) then	
				
					factor=(48.0d0*Epsi/Sigma*2.0d0)*( (Sigma/modraio)**14.0d0-0.5d0*(Sigma/modraio)**8.0d0)
					fx(i)=fx(i)+factor*dr(1)
					fy(i)=fy(i)+factor*dr(2)
					fz(i)=fz(i)+factor*dr(3)

					fx(j)=fx(j)-factor*dr(1) !3a lei
					fy(j)=fy(j)-factor*dr(2)
					fz(j)=fz(j)-factor*dr(3)
				
				end if
				j = j + 1
			end do
			i = i + 1
		end do
		else if (switch == 2) then
		utot = 0.0d0    !zerar variável de acumulação
		i = 1
		
		do while (i.le.n-1)
			j = i + 1
			do while (j.le.n)	
				dr(1) = x(i) - x(j)
				dr(2) = y(i) - y(j)
				dr(3) = z(i) - z(j)
				dr(1) = dr(1) - L * ANINT(dr(1) / L)
				dr(2) = dr(2) - L * ANINT(dr(2) / L)
				dr(3) = dr(3) - L * ANINT(dr(3) / L)
				modraio=SQRT( dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3))
				if (modraio.le.raiodecorte) then	
					if (pot == 2) then
						utot = utot + 4.0d0*Epsi*( (Sigma/modraio)**12.0d0-(Sigma/modraio)**6.0d0)+Epsi  !aqui calcular a energia 					
					else 
    	        		utot = utot + 4.0d0*Epsi*( (Sigma/modraio)**12.0d0-(Sigma/modraio)**6.0d0)
					end if
				end if
				j = j + 1
			end do
			i = i + 1
		end do				
	
		else 
			print*,"deu ruim"
			stop

		end if

	end subroutine force

	subroutine stats (N,X,Xmed,DesvP)
		integer::N  !numero de amostras
		real,dimension(N) :: X !vetor contendo valores de amostras
		real(8) :: Xmed
		real(8) :: DesvP
		integer :: I
		DesvP = 0.0d0
		Xmed = 0.0d0
		do I = 1 , N
			Xmed = Xmed + X(I)
		end do
		Xmed = Xmed/real(N)	
		do I = 1 , N
			DesvP = DesvP + ( X(I) - Xmed)**2
		end do
		DesvP = SQRT(DesvP/(N - 1))
	end subroutine






		
	FUNCTION ran2(idum)
		  INTEGER idum,idum2
		  INTEGER IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
		  real(8) ran2,AM,EPS,RNMX
		  PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,&
		       IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,&
		       NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
		  INTEGER j,k,iv(NTAB),iy
		  SAVE iv,iy,idum2
		  DATA idum2/123456789/, iv/NTAB*0/, iy/0/
		  if (idum.le.0) then
		     idum=max(-idum,1)
		     idum2=idum
		     do j=NTAB+8,1,-1
		        k=idum/IQ1
		        idum=IA1*(idum-k*IQ1)-k*IR1
		        if (idum.lt.0) idum=idum+IM1
		        if (j.le.NTAB) iv(j)=idum
		     enddo
		     iy=iv(1)
		  endif
		  k=idum/IQ1
		  idum=IA1*(idum-k*IQ1)-k*IR1
		  if (idum.lt.0) idum=idum+IM1
		  k=idum2/IQ2
		  idum2=IA2*(idum2-k*IQ2)-k*IR2
		  if (idum2.lt.0) idum2=idum2+IM2
		  j=1+iy/NDIV
		  iy=iv(j)-idum2
		  iv(j)=idum
		  if(iy.lt.1)iy=iy+IMM1
		  ran2=min(AM*iy,RNMX)
		  return
		END FUNCTION ran2





	end program velocity_verlet
