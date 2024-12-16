using Microsoft.AspNetCore.Mvc;
using System.Text.Json;

namespace Homework6;

[ApiController]
[Route("[controller]")]
public class OptionsDataController : ControllerBase
{

    [HttpGet("/OptionsData")]
    public ActionResult<Option> GetOptionsData()
    {
        Console.WriteLine("Getting");
        FinanceContext db = new FinanceContext();

        var all_trades = db.Trades.ToArray();

        List<int> fi_ids = new List<int>();

        foreach (var trade in all_trades)
        {
            // Console.WriteLine(trade.financialinstrumentid);
            fi_ids.Add(trade.financialinstrumentid);
        }

        fi_ids.ForEach(x => Console.WriteLine(x));

        var result = db.EuropeanOptions.Where(x => fi_ids.Contains(x.Id)).ToArray();
        var result2 = db.DigitalOptions.Where(x => fi_ids.Contains(x.Id));
        var result3 = db.AsianOptions.Where(x => fi_ids.Contains(x.Id));
        var result4 = db.RangeOptions.Where(x => fi_ids.Contains(x.Id));
        var result5 = db.LookbackOptions.Where(x => fi_ids.Contains(x.Id));
        var result6 = db.BarrierOptions.Where(x => fi_ids.Contains(x.Id));

        List<Dictionary<string, object>> mylist = new List<Dictionary<string, object>>();

        foreach (var res in result) {
            Dictionary<string, object> toadd = new Dictionary<string, object>();
            toadd["European"] = res;
            mylist.Add(toadd);
        }

        foreach (var res in result2) {
            Dictionary<string, object> toadd = new Dictionary<string, object>();
            toadd["Digital"] = res;
            mylist.Add(toadd);
        }

        foreach (var res in result3) {
            Dictionary<string, object> toadd = new Dictionary<string, object>();
            toadd["Asian"] = res;
            mylist.Add(toadd);
        }

        foreach (var res in result4) {
            Dictionary<string, object> toadd = new Dictionary<string, object>();
            toadd["Range"] = res;
            mylist.Add(toadd);
        }

        foreach (var res in result5) {
            Dictionary<string, object> toadd = new Dictionary<string, object>();
            toadd["Lookback"] = res;
            mylist.Add(toadd);
        }

        foreach (var res in result6) {
            Dictionary<string, object> toadd = new Dictionary<string, object>();
            toadd["Barrier"] = res;
            mylist.Add(toadd);
        }

        return Ok(mylist.ToArray());
    }
}
