using Microsoft.AspNetCore.Mvc;
using System.Text.Json;

namespace Homework6;

[ApiController]
[Route("[controller]")]
public class UnitsController : ControllerBase
{

    [HttpPost("/Units")]
    public ActionResult<Unit> PostUnits([FromBody] Unit units)
    {
        Console.WriteLine("Post");
        modifydb.modifyunits(units.Id, units.typeUnit, units.sizeUnit);
        return Ok(units);
    }

    [HttpGet("/Units")]
    public ActionResult<Unit> GetUnits()
    {
        Console.WriteLine("uning get units");
        FinanceContext db = new FinanceContext();
        Console.WriteLine(db);
        return Ok(db.Units.ToArray());
    }
}
