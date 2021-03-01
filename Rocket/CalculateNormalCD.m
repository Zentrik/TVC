function NormalCD = CalculateNormalCD(K, AoA, Reference_Area, NoseCone_PlanformArea, BodyTube_PlanformArea)
    NormalCD = sin(AoA) * ( sin(AoA) * K * (NoseCone_PlanformArea + BodyTube_PlanformArea) / Reference_Area + 2);