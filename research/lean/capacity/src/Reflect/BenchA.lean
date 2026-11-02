import Reflect.Core0
namespace HodgeComputer.Reflect
set_option maxRecDepth 100000
-- generation only (two 100-term polys)
example : (mkPoly 1 100).length + (mkPoly 2 100).length = 200 := by decide +kernel
end HodgeComputer.Reflect
