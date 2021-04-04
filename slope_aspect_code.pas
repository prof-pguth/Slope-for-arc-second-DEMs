This is extracted from the Delphi (Object Pascal) source code for MICRODEM, which totals around 400,000 lines of code


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


         procedure GetAspect( var SlopeAsp : tSlopeAspectRec);   inline;
         begin
            if (abs(SlopeAsp.dzdx) < 0.001) and (abs(SlopeAsp.dzdy) < 0.001) then begin
               SlopeAsp.Dir := cdFlat;
               SlopeAsp.AspectDir := MaxSmallInt;
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

               if (SlopeAsp.z < SlopeAsp.zne) and (SlopeAsp.z < SlopeAsp.znw) and (SlopeAsp.z < SlopeAsp.zn) and (SlopeAsp.z < SlopeAsp.ze) and
                  (SlopeAsp.z < SlopeAsp.zw) and (SlopeAsp.z < SlopeAsp.zse) and (SlopeAsp.z < SlopeAsp.zsw) and (SlopeAsp.z < SlopeAsp.zs) then begin
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
   Lat,Long : float64;
begin
   Result := false;
   if SurroundedPointElevs(Col,Row,SlopeAsp.znw,SlopeAsp.zw,SlopeAsp.zsw,SlopeAsp.zn,SlopeAsp.z,SlopeAsp.zs,SlopeAsp.zne,SlopeAsp.ze,SlopeAsp.zse,MDDef.SlopeRegionRadius) then begin
      Result := true;
      PixelSpacingAndRotation(Col,Row,Lat,Long,SlopeAsp.dx,SlopeAsp.dy,SlopeAsp.GridTrueAngle);
      SlopeAsp.dzdx := 1 / SlopeAsp.dx / 6 * (+SlopeAsp.zne+SlopeAsp.ze+SlopeAsp.zse-SlopeAsp.znw-SlopeAsp.zw-SlopeAsp.zsw);
      SlopeAsp.dzdy := 1 / SlopeAsp.dy / 6 * (+SlopeAsp.znw+SlopeAsp.zn+SlopeAsp.zne-SlopeAsp.zsw-SlopeAsp.zs-SlopeAsp.zse);
      SlopeAsp.Slope := sqrt(sqr(SlopeAsp.dzdx) + sqr(SlopeAsp.dzdy));
      SlopeAsp.SlopePercent := 100 * SlopeAsp.Slope;
      SlopeAsp.SlopeDegree := ArcTan(SlopeAsp.Slope) / DegToRad;
      GetAspect(SlopeAsp);
   end
   else begin
      SlopeAsp.Slope := MaxSmallInt;
      SlopeAsp.AspectDir := MaxSmallInt;
      SlopeAsp.Slope := 0;
   end;
end {proc};


function HeadingOfLine(dx,dy : float64) : float64;
{returns the heading of a line, with the value between 0 and 360 in degrees}
{   values start with 0 north and proceed clockwise                        }
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


procedure tDEMDataSet.PixelSpacingAndRotation(Col,Row : integer; var Lat,Long,xdistance,ydistance,GridTrueAngle : float64);
var
   Lat1,Long1,Lat2,Long2,Bearing : float64;
begin
//this calls routines to deal with all possible cases of the DEM geometry
   DEMGridToLatLongDegree(Col,Row,Lat,Long);
   DEMGridToLatLongDegree(Col,pred(Row),Lat1,Long1);
   DEMGridToLatLongDegree(Col,succ(Row),Lat2,Long2);
   VincentyCalculateDistanceBearing(Lat1,Long1,Lat2,Long2,yDistance,GridTrueAngle);
   if GridTrueAngle > 180 then GridTrueAngle := - (360 - GridTrueAngle);

   DEMGridToLatLongDegree(pred(Col),Row,Lat1,Long1);
   DEMGridToLatLongDegree(succ(Col),Row,Lat2,Long2);
   VincentyCalculateDistanceBearing(Lat1,Long1,Lat2,Long2,xDistance,Bearing);
   xdistance := 0.5 * xDistance;
   ydistance := 0.5 * yDistance;
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


