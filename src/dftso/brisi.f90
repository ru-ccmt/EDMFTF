program brisi
  COMPLEX*16 :: a(3), b(3)

  a(1) = dcmplx(1,1)
  a(2) = dcmplx(0,1)
  a(3) = dcmplx(1,2)
  b(1) = dcmplx(3,1)
  b(2) = dcmplx(1,7)
  b(3) = dcmplx(0,6)

  print *, dot_product(a,b)
  print *, zdotc(3,a,1,b,1)
  
end program brisi
