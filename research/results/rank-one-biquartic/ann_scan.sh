for T in 'X1**3' 'X2**3' '(X0+X1)**3' '(X0+2*X1-X2)**3' 'Y1**3' 'Y2**3' '(Y0+Y1)**3' 'X0**2' 'X0' 'X0**4' 'X0**2*X1' 'X0*Y0' 'X0**2*Y0' 'X0**3*Y0' 'X0**3*Y0**2' 'L' 'M'; do
  if [ "$T" = "L" ]; then T='X0**2 + X0*X1 + 2*X1*X2 + 3*X2**2'; fi
  if [ "$T" = "M" ]; then T='Y0**2 + Y0*Y1 + 2*Y1*Y2 + 3*Y2**2'; fi
  echo "== TARGET $T"
  TARGET="$T" python3 e1kernel_gen.py 5 6 7 8 2>/dev/null | grep -o '"d": [0-9]*, "parD": [0-9]*.*"kernel_mod_J": [0-9]*' | sed 's/"rows.*"dim_E1_space"/"E1"/'
done
