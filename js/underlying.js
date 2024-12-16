document.onload = handle_get();

function handle_submit() {

    var request = new XMLHttpRequest();
    request.open("POST", "http://localhost:5115/Underlying");
    request.setRequestHeader("Content-Type", "application/json");

    request.onreadystatechange = function () {
        if (request.readyState === 4 && request.status === 200) {
            var output = JSON.parse(request.responseText);
            console.log(output);
        }
    };

    var object = {};
    var form_data = new FormData(document.forms.form_underlying);
    console.log(form_data);

    form_data.forEach(function(value, key) {
        if (value === 'True') {
            value = true;
        } 

        if (value === 'False') {
            value = false;
        }

        object[key] = value;
    });

    var json = JSON.stringify(object);

    console.log(json)
    request.send(json);
}

const button = document.getElementById("underlying_submit");

button.addEventListener("click", handle_submit);

function handle_get() {

    var request = new XMLHttpRequest();
    request.open("Get", "http://localhost:5115/Underlying");
    request.setRequestHeader("Content-Type", "application/json");

    request.onreadystatechange = function () {
        if (request.readyState === 4 && request.status === 200) {
            var output = JSON.parse(request.responseText);
            console.log(output);
            var mytable = document.getElementById("underlyingtable");
            let rowcount = mytable.rows.length;
            for (let i = 1; i < rowcount; i++) {
                mytable.deleteRow(-1);
            }
            for (let i = 0; i < output.length; i++) {
                var row = mytable.insertRow(mytable.rows.length);
                var idCol = row.insertCell();
                idCol.innerHTML = output[i].id;
                var mymarketid = row.insertCell();
                mymarketid.innerHTML = output[i].marketid;
                var myname = row.insertCell();
                myname.innerHTML = output[i].name;
                var mysymbol = row.insertCell();
                mysymbol.innerHTML = output[i].symbol;
            }
        }
    };
    console.log("JS Get");
    request.send();
}

const button_get = document.getElementById("underlying_refresh");
button_get.addEventListener("click", handle_get);