# Attitude

The attitude module contains function for common attitude operations, different
attitude representations (Quaternions, Euler Angles, Euler Axis & Angle, Rotation 
Matrices), and transformations between different attitude representations.

```@docs
Rx
Ry
Rz
```

## Quaternions

```@docs
Quaternion
as_vector(q::Quaternion)
as_matrix(q::Quaternion)
copy(q::Quaternion)
deepcopy(q::Quaternion)
slerp
```

## EulerAngle

```@docs
EulerAngle
as_vector(e::EulerAngle)
as_matrix(e::EulerAngle)
copy(e::EulerAngle)
deepcopy(e::EulerAngle)
```

## EulerAxis

```@docs
EulerAxis
as_vector(e::EulerAxis)
as_matrix(e::EulerAxis)
copy(e::EulerAxis)
deepcopy(e::EulerAxis)
```