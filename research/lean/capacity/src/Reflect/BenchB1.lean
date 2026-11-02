import Reflect.Core1
namespace HodgeComputer.Reflect1
set_option maxRecDepth 100000
example : 0 < (mul (mkPoly 1 100) (mkPoly 2 100)).length := by decide +kernel
end HodgeComputer.Reflect1
