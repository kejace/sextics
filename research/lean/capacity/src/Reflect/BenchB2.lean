import Reflect.Core2
namespace HodgeComputer.Reflect2
set_option maxRecDepth 100000
example : 0 < (mul (mkPoly 1 100) (mkPoly 2 100)).length := by decide +kernel
end HodgeComputer.Reflect2
