# Taste Fortran

* Quickstart tutorial - https://fortran-lang.org/en/learn/quickstart
* Fortran90 tutorial - https://web.stanford.edu/class/me200c/tutorial_90
* Fortran Tutorial - http://seismic.yonsei.ac.kr/fortran/index.html
* fortran sample codes - https://people.math.sc.edu/Burkardt/f_src/f_src.html
* f77 to f90 - https://people.math.sc.edu/Burkardt/f_src/f77_to_f90/f77_to_f90.html
* 강좌
    * https://shlee1990.tistory.com/category/%ED%94%84%EB%A1%9C%EA%B7%B8%EB%9E%98%EB%B0%8D%20%EC%96%B8%EC%96%B4/Fortran
    * https://bd.kma.go.kr/kma2020/dta/edu/KBP57200_Fortran.do?pageNum=5&menuCd=F040304000
* etc
    * 변수선언 - https://pcium.tistory.com/12?category=741185
    * http://www.personal.psu.edu/hdk/fortran.html
    * https://jblevins.org/mirror/amiller

## Compile
```powershll
# gfortran -fallow-argument-mismatch -std=legacy -o out.exe hello.f
# gfortran -fallow-argument-mismatch -std=legacy -o out.exe hello.f90
gfortran -std=legacy -o out.exe hello.f
gfortran -std=legacy -o out.exe hello.f90
```

## Convert f77 to f90
```powershell
gfortran -o f77tof90.exe .\f77_to_f90.f90
```
