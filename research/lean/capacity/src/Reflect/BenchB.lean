import Reflect.Core0
namespace HodgeComputer.Reflect
set_option maxRecDepth 100000
-- one product 100x100, check only its length is positive
example : 0 < (mul (mkPoly 1 100) (mkPoly 2 100)).length := by decide +kernel
end HodgeComputer.Reflect
