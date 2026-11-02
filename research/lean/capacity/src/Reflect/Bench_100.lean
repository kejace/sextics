import Reflect.Core0
namespace HodgeComputer.Reflect
set_option maxRecDepth 100000
-- commutativity check forces two full products of 100 x 100 terms
example : mul (mkPoly 1 100) (mkPoly 2 100) = mul (mkPoly 2 100) (mkPoly 1 100) := by decide +kernel
end HodgeComputer.Reflect
