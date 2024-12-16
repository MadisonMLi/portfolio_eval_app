document.onload = handle_get();

function handle_submit() {

    var request = new XMLHttpRequest();
    request.open("POST", "http://localhost:5115/Units");
    request.setRequestHeader("Content-Type", "application/json");

    request.onreadystatechange = function () {
        if (request.readyState === 4 && request.status === 200) {
            var output = JSON.parse(request.responseText);
            console.log(output);
        }
    };

    var object = {};
    var form_data = new FormData(document.forms.form_units);
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

const button = document.getElementById("units_submit");

button.addEventListener("click", handle_submit);

function handle_get() {

    var request = new XMLHttpRequest();
    request.open("Get", "http://localhost:5115/Units");
    request.setRequestHeader("Content-Type", "application/json");

    request.onreadystatechange = function () {
        if (request.readyState === 4 && request.status === 200) {
            var output = JSON.parse(request.responseText);
            console.log(output);
            var mytable = document.getElementById("unitstable");
            let rowcount = mytable.rows.length;
            for (let i = 1; i < rowcount; i++) {
                mytable.deleteRow(-1);
            }
            for (let i = 0; i < output.length; i++) {
                var row = mytable.insertRow(mytable.rows.length);
                var idCol = row.insertCell();
                idCol.innerHTML = output[i].id;
                var mytype = row.insertCell();
                mytype.innerHTML = output[i].typeUnit;
                var mysize = row.insertCell();
                mysize.innerHTML = output[i].sizeUnit;
            }
        }
    };

    console.log("gg");
    request.send();

}

const button_get = document.getElementById("units_refresh");
button_get.addEventListener("click", handle_get);