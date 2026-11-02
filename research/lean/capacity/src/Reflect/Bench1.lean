import Reflect.Core0
namespace HodgeComputer.Reflect
set_option maxRecDepth 100000
-- correctness smoke: (p*q) computed two ways must agree: mul p q = mul q p
example : mul (mkPoly 1 30) (mkPoly 2 30) = mul (mkPoly 2 30) (mkPoly 1 30) := by decide +kernel
end HodgeComputer.Reflect
