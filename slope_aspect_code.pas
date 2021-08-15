//This is extracted from the Delphi (Object Pascal) source code for MICRODEM, which totals around 400,000 lines of code
//If you are serious about working with this in Delphi, contact me, pguth@usna.edu or pguth@verizon.net
//this version is from 15 August 2021


type
   tSlopeAspectRec = record
       z,znw,zw,zsw,zn,zs,zne,ze,zse,
       dzdx,dzdy,
       dx,dy,GridTrueAngle,
       Slope,
       SlopePercent,
       SlopeDegree,
       AspectDir   : float64;
       Dir : tCompassDirection;
   end;



function tDEMDataSet.GetSlopeAndAspect(Col,Row : integer; var SlopeAsp : tSlopeAspectRec) : boolean;
//based on coordinates in the DEM grid, which starts in the SW corner
type
   ShortReal = array[1..4] of float64;
var
   sl     : array[1..4] of float64;
   AspDir : array[2..4] of tCompassDirection;
   UseDiaSpace : float64;


         procedure AddSlope(z,z1,z2,Dist : float64);
         var
            d1,d2 : float64;
         begin
            d1 := (z - z1);
            d2 := (z - z2);
            d1 := abs(d1);
            SlopeAsp.Slope := SlopeAsp.Slope + d1 / Dist;
            d2 := abs(d2);
            SlopeAsp.Slope := SlopeAsp.Slope + d2 / Dist;
         end;


         procedure MaxDownHill(z,z1,z2,Dist : float64; var Max : float64; Dir1,Dir2 : tCompassDirection; var AspDir : tCompassDirection);    inline;
         var
            d1,d2 : float64;
         begin {proc MaxDownHill}
            d1 := (z - z1);
            d2 := (z - z2);
            if (d1 = d2) then if Odd(Random(100)) then d1 := 0 else d2 := 0;
            if (d1 > d2) and (d1 > 0) then begin
               AspDir := Dir1;
               Max := d1 / Dist;
            end
            else if (d2 > d1) and (d2 > 0) then begin
               AspDir := Dir2;
               Max := d2 / Dist;
            end
            else begin
               Max := 0;
               AspDir := cdPit;
            end;
         end {proc MaxSlope};


         procedure MaxSlope(z,z1,z2,Dist : float64; var Max : float64; Dir1,Dir2 : tCompassDirection; var AspDir : tCompassDirection); inline;
         const
            Opp : array[tCompassDirection] of tCompassDirection = (cdS,cdSW,cdW,cdNW,cdN,cdNE,cdE,cdSE,cdFlat,cdPit);
         var
            d1,d2 : float64;
         begin {proc MaxSlope}
            d1 := (z - z1);
            d2 := (z - z2);
            if (d1 = d2) then if Odd(Random(100)) then d1 := 0 else d2 := 0;

            if abs(d1) > abs(d2) then begin
               if d1 < 0 then AspDir := Opp[Dir1] else AspDir := Dir1;
               d1 := abs(d1);
               Max := d1 / Dist;
            end
            else begin
               if d2 < 0 then AspDir := Opp[Dir2] else AspDir := Dir2;
               d2 := abs(d2);
               Max := d2 / Dist;
            end;
         end {proc MaxSlope};


         procedure MaxSlopeComputations;
         begin
             MaxSlope(SlopeAsp.z,SlopeAsp.zs,SlopeAsp.zn,SlopeAsp.dy,sl[1],cdS,cdN,SlopeAsp.Dir);
             MaxSlope(SlopeAsp.z,SlopeAsp.zw,SlopeAsp.ze,SlopeAsp.dx,sl[2],cdW,cdE,AspDir[2]);
             MaxSlope(SlopeAsp.z,SlopeAsp.zsw,SlopeAsp.znw,UseDiaSpace,sl[3],cdSW,cdNW,AspDir[3]);
             MaxSlope(SlopeAsp.z,SlopeAsp.zse,SlopeAsp.zne,UseDiaSpace,sl[4],cdSE,cdNE,AspDir[4]);
         end;


         function IsPit(var SlopeAsp : tSlopeAspectRec) : boolean; inline;
         begin
            Result := (SlopeAsp.z < SlopeAsp.zne) and (SlopeAsp.z < SlopeAsp.znw) and (SlopeAsp.z < SlopeAsp.zn) and (SlopeAsp.z < SlopeAsp.ze) and
                  (SlopeAsp.z < SlopeAsp.zw) and (SlopeAsp.z < SlopeAsp.zse) and (SlopeAsp.z < SlopeAsp.zsw) and (SlopeAsp.z < SlopeAsp.zs);
         end;


         procedure GetAspect(var SlopeAsp : tSlopeAspectRec); inline;
         begin
            if (abs(SlopeAsp.dzdx) < 0.001) and (abs(SlopeAsp.dzdy) < 0.001) then begin
               SlopeAsp.AspectDir := MaxSmallInt;
               SlopeAsp.Dir := cdFlat;
            end
            else begin
               //modified atan2 function
               //   standard math convention puts 0 on the x axis and angles increases counterclockwise,
               //   use geographic conventions, where N (0) is on the y axis, and angles increase clockwise
               SlopeAsp.AspectDir := HeadingOfLine(SlopeAsp.dzdx,SlopeAsp.dzdy) + 180;
               //correct so that aspect is referenced to true north
               SlopeAsp.AspectDir := SlopeAsp.AspectDir + SlopeAsp.GridTrueAngle;
               if (SlopeAsp.AspectDir > 360) then SlopeAsp.AspectDir := SlopeAsp.AspectDir - 360;
               if (SlopeAsp.AspectDir < 0) then SlopeAsp.AspectDir := SlopeAsp.AspectDir + 360;

               if IsPit(SlopeAsp) then begin
                     SlopeAsp.Dir := cdPit;
               end
               else begin
                  if (SlopeAsp.AspectDir < 22.5) or (SlopeAsp.AspectDir > 337.5) then SlopeAsp.Dir := cdN
                  else if (SlopeAsp.AspectDir < 67.5) then SlopeAsp.Dir := cdNE
                  else if (SlopeAsp.AspectDir < 112.5) then SlopeAsp.Dir := cdE
                  else if (SlopeAsp.AspectDir < 157.5) then SlopeAsp.Dir := cdSE
                  else if (SlopeAsp.AspectDir < 202.5) then SlopeAsp.Dir := cdS
                  else if (SlopeAsp.AspectDir < 247.5) then SlopeAsp.Dir := cdSW
                  else if (SlopeAsp.AspectDir < 292.5) then SlopeAsp.Dir := cdW
                  else SlopeAsp.Dir := cdNW;
               end;
            end;
         end;



var
   j : integer;
   Lat,Long : float64;
begin
   Result := false;
   SlopeAsp.AspectDir := MaxSmallInt;
   SlopeAsp.Slope := 0;
   if SurroundedPointElevs(Col,Row,SlopeAsp.znw,SlopeAsp.zw,SlopeAsp.zsw,SlopeAsp.zn,SlopeAsp.z,SlopeAsp.zs,SlopeAsp.zne,SlopeAsp.ze,SlopeAsp.zse,MDDef.SlopeRegionRadius) then with SlopeAsp do begin
      Result := true;

      PixelSpacingAndRotation(Col,Row,Lat,Long,SlopeAsp.dx,SlopeAsp.dy,SlopeAsp.GridTrueAngle,MDDef.QuickSlopeSpacings);
      SlopeAsp.dx := SlopeAsp.dx * MDDef.SlopeRegionRadius;
      SlopeAsp.dy := SlopeAsp.dy * MDDef.SlopeRegionRadius;

      if (MDDef.SlopeAlg in [smEightNeighborsUnweighted,smFourNeighbors,smEightNeighborsWeighted,smEightNeighborsWeightedByDistance,smFrameFiniteDifference,smSimpleDifference,smONeillAndMark]) then begin
         if (MDDef.SlopeAlg in [smEightNeighborsUnweighted]) then begin
            SlopeAsp.dzdx := 1 / SlopeAsp.dx / 6 * (+SlopeAsp.zne+SlopeAsp.ze+SlopeAsp.zse-SlopeAsp.znw-SlopeAsp.zw-SlopeAsp.zsw);
            SlopeAsp.dzdy := 1 / SlopeAsp.dy / 6 * (+SlopeAsp.znw+SlopeAsp.zn+SlopeAsp.zne-SlopeAsp.zsw-SlopeAsp.zs-SlopeAsp.zse);
         end
         else if (MDdef.SlopeAlg = smEightNeighborsWeighted) then begin  //Horn method
            SlopeAsp.dzdx := 0.125 * ( SlopeAsp.zne + (2 * SlopeAsp.ze) + SlopeAsp.zse  -SlopeAsp.znw - (2 * SlopeAsp.zw) - SlopeAsp.zsw) / SlopeAsp.dx;
            SlopeAsp.dzdy := 0.125 * ( SlopeAsp.znw + (2 * SlopeAsp.zn) + SlopeAsp.zne  -SlopeAsp.zsw - (2 * SlopeAsp.zs) - SlopeAsp.zse) / SlopeAsp.dy;
         end
         else if (MDdef.SlopeAlg = smEightNeighborsWeightedByDistance) then begin
            SlopeAsp.dzdy := 1 / (4 + 2 * sqrt_2) * ((SlopeAsp.znw + Sqrt_2 * SlopeAsp.zn + SlopeAsp.zne) - (SlopeAsp.zsw + Sqrt_2 * SlopeAsp.zs + SlopeAsp.zse)) / SlopeAsp.dy;
            SlopeAsp.dzdx := 1 / (4 + 2 * sqrt_2) * ((SlopeAsp.zne + Sqrt_2 * SlopeAsp.ze + SlopeAsp.zse) - (SlopeAsp.znw + Sqrt_2 * SlopeAsp.zw + SlopeAsp.zsw)) / SlopeAsp.dx;
         end
         else if (MDDef.SlopeAlg = smFrameFiniteDifference) then begin
            SlopeAsp.dzdy := (SlopeAsp.znw - SlopeAsp.zsw + SlopeAsp.zne - SlopeAsp.zse) * 0.25 / SlopeAsp.dy;
            SlopeAsp.dzdx := (SlopeAsp.zse - SlopeAsp.zsw + SlopeAsp.zne - SlopeAsp.znw) * 0.25 / SlopeAsp.dx;
         end
         else if (MDDef.SlopeAlg = smFourNeighbors) then begin
            SlopeAsp.dzdy := (SlopeAsp.zn - SlopeAsp.zs) * 0.5 / SlopeAsp.dy;
            SlopeAsp.dzdx := (SlopeAsp.ze - SlopeAsp.zw) * 0.5 / SlopeAsp.dx;
         end
         else if (MDDef.SlopeAlg = smSimpleDifference) then begin
            SlopeAsp.dzdy := (SlopeAsp.z - SlopeAsp.zs) * 0.5 / SlopeAsp.dy;
            SlopeAsp.dzdx := (SlopeAsp.z - SlopeAsp.zw) * 0.5 / SlopeAsp.dx;
         end
         else if (MDDef.SlopeAlg in [smONeillAndMark]) then begin
            SlopeAsp.dzdy := (SlopeAsp.zn - SlopeAsp.z) / SlopeAsp.dy;
            SlopeAsp.dzdx := (SlopeAsp.ze - SlopeAsp.z) / SlopeAsp.dx;
         end;
         SlopeAsp.Slope := sqrt(sqr(SlopeAsp.dzdx) + sqr(SlopeAsp.dzdy));
         SlopeAsp.SlopePercent := 100 * SlopeAsp.Slope;
         SlopeAsp.SlopeDegree := ArcTan(SlopeAsp.Slope) / DegToRad;
         GetAspect(SlopeAsp);
      end
      else begin
         UseDiaSpace := DiagSpaceByDEMrow^[Row] * MDDef.SlopeRegionRadius;
         if (MDDef.SlopeAlg in [smGuthHybrid]) then begin  //no longer recommended
             MaxSlopeComputations;
             for j := 2 to 4 do if sl[j] > sl[1] then sl[1] := sl[j];
             Slope := Sl[1];
             dzdx := (SlopeAsp.zne + SlopeAsp.ze + SlopeAsp.zse - SlopeAsp.zsw - SlopeAsp.zw - SlopeAsp.znw) / 6 / SlopeAsp.dx;
             dzdy := (SlopeAsp.znw + SlopeAsp.zn + SlopeAsp.zne - SlopeAsp.zsw - SlopeAsp.zs - SlopeAsp.zse) / 6 / SlopeAsp.dy;
             GetAspect(SlopeAsp);
         end
         else if (MDDef.SlopeAlg in [smSteepestNeighbor,smAverageNeighbor]) then begin
            MaxSlopeComputations;
            for j := 2 to 4 do if sl[j] > sl[1] then begin
                sl[1] := sl[j];
                Dir := AspDir[j];
            end {for j};
            Slope := Sl[1];
            AspectDir := 45.0 * Ord(Dir);
            if (MDDef.SlopeAlg = smAverageNeighbor) then begin
               Slope := 0;
               AddSlope(z,zs,zn,SlopeAsp.dy);
               AddSlope(z,zw,ze,SlopeAsp.dx);
               AddSlope(z,zsw,znw,UseDiaSpace);
               AddSlope(z,zse,zne,UseDiaSpace);
               Slope := Slope / 8;
            end;
         end
         else if (MDDef.SlopeAlg = smMaxDownhillSlope) then begin
            MaxDownHill(z,zs,zn,SlopeAsp.dy,sl[1],cdS,cdN,Dir);
            MaxDownHill(z,zw,ze,SlopeAsp.dx,sl[2],cdW,cdE,AspDir[2]);
            MaxDownHill(z,zsw,znw,UseDiaSpace,sl[3],cdSW,cdNW,AspDir[3]);
            MaxDownHill(z,zse,zne,UseDiaSpace,sl[4],cdSE,cdNE,AspDir[4]);
            for j := 2 to 4 do if sl[j] > sl[1] then begin
                sl[1] := sl[j];
                Dir := AspDir[j];
            end {for j};
            Slope := Sl[1];
            AspectDir := 45.0 * Ord(Dir);
         end;
         if IsPit(SlopeAsp) then SlopeAsp.Dir := cdPit;
      end;
      SlopeAsp.SlopePercent := 100 * Slope;
      SlopeAsp.SlopeDegree := ArcTan(Slope) / DegToRad;
   end;
end {proc GetSlopeAndAspect};



function HeadingOfLine(dx,dy : float64) : float64; inline;
{returns the heading of a line, with the value between 0 and 360 in degrees}
{   values start with 0 north and proceed clockwise, geographic convention }
{modified atan2 function, deals with 0 values in dx or dy                  }
begin
   if abs(dy) < 0.000001 then
      if (dx > 0) then Result := 90 else Result := 270
   else
      if (abs(dx) < 0.00001) then
         if dy < 0 then Result := 180 else Result := 0
      else begin
         Result := ArcTan(dx / dy) / DegToRad;
         if (dx < 0) and (dy < 0) then Result := Result + 180;
         if (dx > 0) and (dy < 0) then Result := Result - 180;
         if Result < 0 then Result := Result + 360;
      end;
end {proc HeadingOfLine};



procedure tDEMDataSet.PixelSpacingAndRotation(Col,Row : integer; var Lat,Long,xdistance,ydistance,GridTrueAngle : float64; Quick : boolean = true);
//quick mode uses stored values from the DEM initial opening, averages for the DEM which works if the DEM covers a small area
//pixel sizes will agree within about 1 mm for a 1" DEM
//quick mode is about 3 times faster
var
   Lat1,Long1,Lat2,Long2,Lat3,Long3,Lat4,Long4,Bearing : float64;
begin
   DEMGridToLatLongDegree(Col,Row,Lat,Long);
   if Quick and (DEMheader.DEMUsed = ArcSecDEM) then begin
      xdistance := XSpaceByDEMrow^[Row];
      ydistance := AverageYSpace;
      GridTrueAngle := 0;
   end
   else if Quick and (DEMheader.DEMUsed = UTMbasedDEM) then begin
      xdistance := AverageXSpace;
      ydistance := AverageYSpace;
      GridTrueAngle := AverageGridTrue; //computed at the center of DEM during opening;  this will have problems if the DEM covers a large area near edge of UTM zone
   end
   else begin
      //this calls routines to deal with all possible cases of the DEM geometry for other more unusual projections
      DEMGridToLatLongDegree(Col,pred(Row),Lat1,Long1);
      DEMGridToLatLongDegree(Col,succ(Row),Lat2,Long2);
      VincentyCalculateDistanceBearing(Lat1,Long1,Lat2,Long2,yDistance,GridTrueAngle);
      ydistance := 0.5 * yDistance;
      if (GridTrueAngle > 180) then GridTrueAngle := - (360 - GridTrueAngle);

      DEMGridToLatLongDegree(pred(Col),Row,Lat3,Long3);
      DEMGridToLatLongDegree(succ(Col),Row,Lat4,Long4);
      VincentyCalculateDistanceBearing(Lat3,Long3,Lat4,Long4,xDistance,Bearing);
      xdistance := 0.5 * xDistance;
   end;
end;




procedure VincentyFullCalculateDistanceBearing(StartLat,StartLong,EndLat,EndLong : float64; EllipsoidConstants : tMapProjection;  var Distance,Bearing : float64);
//from T. Vincenty, Survey Review, 23, No 176, p 88-93,1975 to calculate lines ranging from a few cm to nearly 20,000 km, with millimetre accuracy
//verified with Excel spreadsheet from Geoscience Australia, http://www.auslig.gov.au/geodesy/datums/vincenty.xls  (http://www.auslig.gov.au/geodesy/datums/calcs.htm)
//requires a and f for the ellipsoid, for WGS84
//   a := 6378137;
//   h_f := 298.2572236;
//input in degrees


   {$IfDef RecordGeotdeticCalcDetailed}
         procedure RecordFirstValues;
         begin
            WriteLineToDebugFile('CalculateDistanceBearing');
            WriteLineToDebugFile('lambda1=' + realToString(lambda1,18,8));
            WriteLineToDebugFile('U1=' + realToString(u1,18,8));
            WriteLineToDebugFile('sinU1=' + realToString(sinu1,18,8));
            WriteLineToDebugFile('cosU1=' + realToString(cosu1,18,8));
            WriteLineToDebugFile('U2=' + realToString(u2,18,8));
            WriteLineToDebugFile('sinU2=' + realToString(sinu2,18,8));
            WriteLineToDebugFile('cosU2=' + realToString(cosu2,18,8));
            WriteLineToDebugFile('realf=' + realToString(realf,18,8));
         end;

         procedure RecordSecondValues
            WriteLineToDebugFile('');
            WriteLineToDebugFile('First lambda: ' + RealToString(Lambda1,18,-8));
            WriteLineToDebugFile('');
            WriteLineToDebugFile('sin_sq_sigma1=' + realToString(sin_sq_sigma1,18,8));
            WriteLineToDebugFile('cos_sigma1=' + realToString(cos_sigma1,18,8));
            WriteLineToDebugFile('tan_sigma1=' + realToString(tan_sigma1,18,8));
            WriteLineToDebugFile('sigma1=' + realToString(Sigma1,18,8));
            WriteLineToDebugFile('two_sigma_m1=' + realToString(two_sigma_m1,18,8));
            WriteLineToDebugFile('cos_two_sigma_m1=' + realToString(cos_two_sigma_m1,18,8));
            WriteLineToDebugFile('Sin_Alpha1=' + realToString(Sin_Alpha1,18,8));
            WriteLineToDebugFile('U_sq1=' + realToString(u_sq1,18,8));
            WriteLineToDebugFile('A1=' + realToString(A1,18,8));
            WriteLineToDebugFile('B1=' + realToString(B1,18,8));
            WriteLineToDebugFile('delta_sigma1=' + realToString(delta_sigma1,18,8));
            WriteLineToDebugFile('c1=' + realToString(c1,18,8));
         end;


         procedure RecordThirdValues;
         begin
            WriteLineToDebugFile('');
            WriteLineToDebugFile('sin_sq_sigma=' + realToString(sin_sq_sigma,18,12));
            WriteLineToDebugFile('cos_sigma=' + realToString(cos_sigma,18,12));
            WriteLineToDebugFile('tan_sigma=' + realToString(tan_sigma,18,12));
            WriteLineToDebugFile('sin_sigma=' + realToString(sin_sigma,18,12));
            WriteLineToDebugFile('sigma=' + realToString(Sigma,18,12));
            WriteLineToDebugFile('Sin_Alpha=' + realToString(Sin_Alpha,18,12));
            WriteLineToDebugFile('cos_two_sigma_m=' + realToString(cos_two_sigma_m,18,12));
            WriteLineToDebugFile('two_sigma_m=' + realToString(two_sigma_m,18,12));
            WriteLineToDebugFile('U_sq=' + realToString(u_sq,18,12));
            WriteLineToDebugFile('A=' + realToString(A,18,12));
            WriteLineToDebugFile('B=' + realToString(B,18,12));
            WriteLineToDebugFile('B_ell=' + realToString(B,18,12));
            WriteLineToDebugFile('delta_sigma=' + realToString(delta_sigma,18,12));
            WriteLineToDebugFile('c=' + realToString(c,18,12));
            WriteLineToDebugFile('lambda=' + realToString(lambda,18,12));
            WriteLineToDebugFile('delta_lambda=' + realToString(abs(Lambda - LastLambda),18,12));
            WriteLineToDebugFile('');
         end;

         procedure RecordFinalValues;
         begin
            WriteLineToDebugFile('Final');
            WriteLineToDebugFile('sin(lambda): ' + realToString(sin(Lambda),18,8)) ;
            WriteLineToDebugFile('Part 1: ' + realToString(cosU2 * sin(Lambda) * cosU2 * sin(Lambda),18,8)) ;
            WriteLineToDebugFile('Part 2: ' + realToString( sqr(cosU1 * sinU2 - sinU1 * cosU2 * cos(Lambda)),18,8) );
            WriteLineToDebugFile('sin_sq_sigma=' + realToString(sin_sq_sigma,18,8));
            WriteLineToDebugFile('cos_sigma=' + realToString(cos_sigma,18,8));
            WriteLineToDebugFile('tan_sigma=' + realToString(tan_sigma,18,8));
            WriteLineToDebugFile('sigma=' + realToString(Sigma,18,8));
            WriteLineToDebugFile('two_sigma_m=' + realToString(two_sigma_m,18,8));
            WriteLineToDebugFile('Sin_Alpha=' + realToString(Sin_Alpha,18,8));
            WriteLineToDebugFile('U_sq=' + realToString(u_sq,18,8));
            WriteLineToDebugFile('A=' + realToString(A,18,8));
            WriteLineToDebugFile('B=' + realToString(B,18,8));
            WriteLineToDebugFile('lambda=' + realToString(lambda,18,8));
            WriteLineToDebugFile('c=' + realToString(c,18,8));
            WriteLineToDebugFile('tan_alpha =' + realToString(tan_alpha ,18,8));
            WriteLineToDebugFile('');
         end;
   {$EndIf}


var
   sin_sq_sigma1,sin_sigma1,cos_sigma1,tan_sigma1,sigma1,sin_alpha1,alpha1,cos_two_sigma_m1,two_sigma_m1,u_sq1,
   A1,B1,delta_sigma1,c1,Lambda1,b_ell,u_sq,tanU1,U1,sinU1,cosU1,A,B,sigma,tanU2,U2,sinU2,cosU2,
   One_minus_f,Lambda,LastLambda,c,alpha,RealF,tan_sigma,sin_alpha,tan_alpha,Lat1,Lat2,Long1,Long2,sin_sq_sigma,
   cos_two_sigma_m,two_sigma_m,delta_sigma,cos_sigma,sin_sigma  : extended;
   i : integer;
begin
   {$IfDef RecordGeotdeticCalc}  WriteLineToDebugFile('FullCalculateDistanceBearing in  Start:  ' + LatLongDegreeToString(StartLat,StartLong) + ' End:    ' + LatLongDegreeToString(EndLat,EndLong));  {$EndIf}

   if (abs(StartLat-EndLat) < 0.0000001) and (abs(StartLong-EndLong) < 0.0000001) then begin
      Distance := 0;
      Bearing := 0;
      exit;
   end;

   Lat1 := StartLat * DegToRad;
   Long1 := StartLong * DegToRad;
   Lat2 := EndLat * DegToRad;
   Long2 := EndLong * DegToRad;

   RealF := 1 / EllipsoidConstants.h_f;
   One_Minus_f := 1 - RealF;
   b_ell := EllipsoidConstants.a * One_Minus_f;

   tanU1 := One_Minus_f * Math.tan(Lat1);
   U1 := arctan(tanU1);
   sinU1 := sin(u1);
   cosU1 := cos(U1);

   tanU2 := One_Minus_f * Math.tan(Lat2);
   U2 := arctan(tanU2);
   sinU2 := sin(u2);
   cosU2 := cos(U2);

   lambda1 := Long2 - Long1;

   {$IfDef RecordGeotdeticCalcDetailed} RecordFirstValues; {$EndIf}

   //formula 14
   sin_sq_sigma1 := (cosU2 * sin(Lambda1) * cosU2 * sin(Lambda1)) + sqr(cosU1 * sinU2 - sinU1 * cosU2 * cos(Lambda1));
   sin_sigma1 := sqrt(sin_sq_sigma1);
   //formula 15
   cos_sigma1 := sinU1 * sinU2 + cosU1 * cos(u2) * cos(Lambda1);
   //formula 16
   tan_sigma1 := sin_sigma1 / cos_sigma1;
   sigma1 := Math.ArcTan2(sin_sigma1, cos_sigma1);
   //formula 17
   if sin_sq_sigma1 < 0.00000000001 then sin_alpha1 := 0
   else sin_alpha1 := cosU1 * cosU2 * sin(lambda1) / sin_sigma1;
   alpha1 := Math.ArcSin(sin_alpha1);
   //formula 18
   cos_two_sigma_m1 := cos_sigma1 - 2 * sinU1 * sinU2 / sqr(cos(alpha1));
   two_sigma_m1 := ArcCos(cos_two_sigma_m1);

   u_sq1 := sqr(cos(Alpha1)) * (sqr(EllipsoidConstants.a) - sqr(b_ell)) / sqr(b_ell);

   //formula 3
   A1 := 1 + u_sq1 / 16384 * (4096 + u_sq1 * ( -768 + u_sq1 * ( 320 - 175 * u_sq1)));
   //formula 4
   B1 := u_sq1 / 1024 * ( 256 + u_sq1 * (-128 + u_sq1 * (74 - 47 * u_sq1)));
   //formula 6
   delta_sigma1 := B1 * sin(sigma1) * (cos(two_sigma_m1) + 0.25 * B1 * ( cos(sigma1) * (-1 + 2 * sqr(cos(two_sigma_m1))
               - 1 / 6 * B1 * cos(two_sigma_m1) * (-3 + 4 * sqr(sqr(sin(sigma1)))) * (-3 + 4 * sqr(sqr(cos(two_sigma_m1)))) )));
   //formula 10
   c1 := RealF  / 16 * sqr(cos(alpha1)) * (4 + RealF * (4 - 3 * sqr(cos(alpha1))));

   //formula 11
   Lambda1 := Lambda1 + (1 - c1) * RealF * sin_alpha1  * (sigma1 + C1 * sin_sigma1 * (cos_two_sigma_m1 + c1 * cos_sigma1 * (-1 + 2 * sqr(cos_two_sigma_m1))));

   {$IfDef RecordGeotdeticCalcDetailed}  RecordSecondValues; {$EndIf}

    Lambda := Lambda1;

   for i := 1 to 1 do begin
   //repeat
      LastLambda := Lambda;
      //formula 14
      sin_sq_sigma := (cosU2 * sin(Lambda) * cosU2 * sin(Lambda)) + sqr(cosU1 * sinU2 - sinU1 * cosU2 * cos(Lambda));
      sin_sigma := sqrt(sin_sq_sigma);
      //formula 15
      cos_sigma := sinU1 * sinU2 + cosU1 * cos(u2) * cos(Lambda);
      //formula 16
      tan_sigma := sin_sigma / cos_sigma;
      sigma := Math.ArcTan2(sin_sigma, cos_sigma);
      //formula 17
      if sin_sq_sigma < 0.00000000001 then sin_alpha := 0
      else sin_alpha := cosU1 * cosU2 * sin(lambda) / sin_sigma;
      alpha := Math.ArcSin(sin_alpha);
      //formula 18
      cos_two_sigma_m := cos_sigma - 2 * sinU1 * sinU2 / sqr(cos(alpha));
      two_sigma_m := ArcCos(cos_two_sigma_m);

      u_sq := sqr(cos(Alpha)) * (sqr(EllipsoidConstants.a) - sqr(b_ell)) / sqr(b_ell);

      //formula 3
      A := 1 + u_sq / 16384 * (4096 + u_sq * ( -768 + u_sq * ( 320 - 175 * u_sq)));
      //formula 4
      B := u_sq / 1024 * ( 256 + u_sq * (-128 + u_sq * (74 - 47 * u_sq)));
      //formula 6
      delta_sigma := B * sin(sigma) * (cos_two_sigma_m + 0.25 * B * (cos(sigma) * (-1 + 2 * sqr(cos_two_sigma_m)) - 1 / 6 * B * cos_two_sigma_m * (-3 + 4 * sqr(sqr(sin(sigma)))) * (-3 + 4 * sqr(sqr(cos_two_sigma_m)) )));
      //formula 10
      c := RealF  / 16 * sqr(cos(alpha)) * (4 + RealF * (4 - 3 * sqr(cos(alpha))));

      //formula 11
      Lambda := (Long2 - Long1) + (1 - c) * RealF * sin_alpha  * (sigma + C * sin_sigma * (cos_two_sigma_m + c * cos_sigma * (-1 + 2 * sqr(cos_two_sigma_m))));

      {$IfDef RecordGeotdeticCalcDetailed} RecordThirdValues; {$EndIf}
   end;

   //formula 20
   tan_alpha := cosU2 * sin(Lambda) / (cosU1*sinU2 - sinU1 * cosU2 * cos(Lambda));

   Bearing := Math.ArcTan2(cosU2 * sin(Lambda), (cosU1*sinU2 - sinU1 * cosU2 * cos(Lambda))) / DegToRad;
   if (Bearing < 0) then Bearing := Bearing + 360;

   {$IfDef RecordGeotdeticCalcDetailed} WriteLineToDebugFile('Call Formula 19'); {$EndIf}

   //formula 19
   Distance := b_ell * (sigma - delta_sigma) * A;

   {$IfDef RecordGeotdeticCalcDetailed} RecordFinalValues; {$EndIf}
   {$IfDef RecordGeotdeticCalc} WriteLineToDebugFile('CalculateDistanceBearing out, Distance=' + RealToString(Distance,12,3) + '  Bearing=' + RealToString(Bearing,8,3)); {$EndIf}
end;



function tDEMDataSet.ReflectanceValue(x,y : integer) : integer;
{returns reflectance value that ranges from 0 (black) to 255 (white); MaxSmallInt if undefined}
{ after Pelton, Colin, 1987, A computer program for hill-shading digital topographic data sets: Computers & Geosciences, vol.13, no.5, p.545-548.}
//https://desktop.arcgis.com/en/arcmap/10.3/tools/spatial-analyst-toolbox/how-hillshade-works.htm
//Hillshade = 255.0 * ((cos(Zenith_rad) * cos(Slope_rad)) + (sin(Zenith_rad) * sin(Slope_rad) * cos(Azimuth_rad - Aspect_rad)))
//Note that if the calculation of the hillshade value is < 0, the output cell value will be = 0.
var
   sum  : float64;
   i : integer;
   SlopeAsp : tSlopeAspectRec;
begin
   if GetSlopeAndAspect(x,y,SlopeAsp) then begin
      Sum := 0;
      for I := 1 to MDdef.UseRefDirs do begin
         Sum := sum + ValidByteRange(round(255.0 * ( (cosDegSunAltitude * cosDeg(RefVertExag * SlopeAsp.SlopeDegree)) + (sinDegSunAltitude * sinDeg(RefVertExag * SlopeAsp.SlopeDegree) * cosDeg(RefPhi[i] - SlopeAsp.AspectDir)))));
      end;
      Result := round(Sum/MDdef.UseRefDirs);
   end
   else begin
      Result := MaxSmallInt;
   end;
end;




