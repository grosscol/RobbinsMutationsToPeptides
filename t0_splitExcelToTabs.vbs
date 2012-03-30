
' Most of this solution is code from the StackOverflow user, brettdj
' http://stackoverflow.com/users/641067/brettdj

if WScript.Arguments.Count = 0 then
  WScript.Echo "Missing parameters"
  WScript.Quit
else
  WScript.Echo "Input: " & WScript.Arguments.Item(0)
end if

Dim strFilename
Dim objFSO

strFilename = WScript.Arguments.Item(0)
Set objFSO = CreateObject("scripting.filesystemobject")

If objFSO.fileexists(strFilename) Then
  Call Writefile(strFilename)
Else
  WScript.echo("No such file!")
  WScript.Quit
End If


Sub Writefile(ByVal strFilename)
  Dim objExcel
  Dim objWB
  Dim objws
  Dim fname
  Set objExcel = CreateObject("Excel.Application")
  Set objWB = objExcel.Workbooks.Open(strFilename)
  For Each objws In objWB.Sheets
    'Make the full file name
    fname = objWB.Path & "\" & objws.Name & ".txt"
    'Check if text file already exists and delete if it does
    If objFSO.FileExists(fname) Then
      objFSO.DeleteFile(fname)
    End If
    objws.Copy()
    objExcel.ActiveWorkbook.SaveAs fname, 21
    WScript.Echo(fname)
    objExcel.ActiveWorkbook.Close False
  Next
  objWB.Close(False)
  objExcel.Quit()
  Set objExcel = Nothing
  Set objFSO = Nothing

End Sub