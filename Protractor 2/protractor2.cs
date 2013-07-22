using System;
using System.Collections.Generic;
using System.Linq;
using UnityEngine;
using System.Text.RegularExpressions;
using System.Reflection;


public class Protractor2Module : PartModule
{
    Rect mainwindow;
    string 
        target = "",
        version = Assembly.GetExecutingAssembly().GetName().Version.ToString();

    #region GUI functions

    public void drawGUI()
    {
        GUI.skin = HighLogic.Skin;
        if (HighLogic.LoadedSceneIsFlight && !FlightDriver.Pause)
        {
            mainwindow = GUILayout.Window(897, mainwindow, mainGUI, "Protractor 2 v" + version, GUILayout.Width(40), GUILayout.Height(20));
        }

    }

    void mainGUI(int windowID)
    {
        try
        {
            GUI.skin = HighLogic.Skin;
            
            CelestialBody planet;

            target = GUILayout.TextField(target, 20);

            try
            {
                foreach (CelestialBody body in FlightGlobals.Bodies)
                {
                    if (body.name == "target") planet = body;
                }
            }
            catch
            {
            }

            GUILayout.BeginHorizontal();
            GUILayout.Label(target + "Data:");
            GUILayout.EndHorizontal();

            GUI.DragWindow();
        }
        catch
        {
            Debug.Log("Error in mainGUI");
        }
    }

    public override void OnStart(PartModule.StartState state)
    {
 	    base.OnStart(state);
        mainwindow = new Rect(Screen.width / 2, Screen.height / 2, 40, 20);
        if (state != StartState.Editor)
        {
            RenderingManager.AddToPostDrawQueue(3, new Callback(drawGUI));
        }
    }

    #endregion


}
