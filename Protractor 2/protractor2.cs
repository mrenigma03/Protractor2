using System;
using System.Collections.Generic;
using System.Linq;
using UnityEngine;
using System.Text.RegularExpressions;
using System.Reflection;


public class Protractor2Module : PartModule
{
    Rect mainwindow;
    string version = Assembly.GetExecutingAssembly().GetName().Version.ToString();
    double UT = 0;

    #region GUI functions

    public void drawGUI()
    {
        GUI.skin = HighLogic.Skin;
        if (HighLogic.LoadedSceneIsFlight && !FlightDriver.Pause)
        {
            mainwindow = GUILayout.Window(897, mainwindow, mainGUI, "Protractor 2 v" + version, GUILayout.Width(150), GUILayout.Height(200));
        }

    }

    void mainGUI(int windowID)
    {
        CelestialBody target;
        ITargetable t;
        

        try
        {
            t = FlightGlobals.fetch.VesselTarget;
            if (t != null && t is CelestialBody)
            {
                target = (CelestialBody) FlightGlobals.fetch.VesselTarget;

                GUILayout.Label("You have selected " + target.name + " as your target.");
                if (UT == 0 || UT < Planetarium.GetUniversalTime())
                {
                    UT = LambertSolver.NextLaunchWindowUT(vessel.mainBody, target);
                }
                GUILayout.Label((UT - Planetarium.GetUniversalTime()).ToString() + " s to window");
                double UT2 = UT;
                Vector3d dv = LambertSolver.EjectionBurn(ref UT2, vessel, target);
                GUILayout.Label((UT2 - Planetarium.GetUniversalTime()).ToString() + " s to window (adj)");

                if (vessel.situation == Vessel.Situations.ORBITING)
                {
                    bool plot = GUILayout.Button("Plot");
                    if (plot)
                    {
                        
                        ManeuverNode mn = vessel.patchedConicSolver.AddManeuverNode(UT);
                        mn.DeltaV = dv;
                        mn.solver.UpdateFlightPlan();
                        plot = false;

                    }
                }
                else
                {
                    GUILayout.Label("Now get into orbit!");
                }

            }
            else
            {
                GUILayout.Label("No valid target selected");
            }
        }
        catch
        {
            GUILayout.Label("Error selecting target");

        }

        

        GUI.DragWindow();
    }

    public override void OnStart(PartModule.StartState state)
    {
 	    base.OnStart(state);
        mainwindow = new Rect(Screen.width / 2, Screen.height / 2, 150, 200);
        if (state != StartState.Editor)
        {
            RenderingManager.AddToPostDrawQueue(3, new Callback(drawGUI));
        }
    }

    #endregion


}
