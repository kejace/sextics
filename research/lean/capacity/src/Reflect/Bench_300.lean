import Reflect.Core0
namespace HodgeComputer.Reflect
set_option maxRecDepth 100000
-- commutativity check forces two full products of 300 x 300 terms
example : mul (mkPoly 1 300) (mkPoly 2 300) = mul (mkPoly 2 300) (mkPoly 1 300) := by decide +kernel
end HodgeComputer.Reflect
